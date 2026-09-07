#ifndef MITM_ROUTER_WORKERS
#define MITM_ROUTER_WORKERS

#include "common.hpp"

namespace mitm {

/******************************** the free stack ********************************/

/* Put a block onto the free stack; the CAS releases everything the pusher did to it.  
 * Reached from Router_Pop, once the receiver has read a block through. 
 */
inline void Router_node::free_push(u32 blk)
{
	u64 w = free_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[blk].store((u32) w, std::memory_order_relaxed);
		u64 w2 = (((w >> 32) + 1) << 32) | blk;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_relaxed))
			return;
	}
}

/* Put a chain the caller linked through blk_link, `first` down to `last`, onto the free stack in one CAS.  
 * Reached from Router_Init (thread 0) and Router_Progress, when the service's stash spills. 
 */
inline void Router_node::free_push_chain(u32 first, u32 last)
{
	u64 w = free_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[last].store((u32) w, std::memory_order_relaxed);
		u64 w2 = (((w >> 32) + 1) << 32) | first;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_relaxed))
			return;
	}
}

/*
 * Pop up to k blocks off the free stack in one CAS, after walking k links.  The generation in the high word
 * changes at every push and pop, so a CAS on a stale view fails: a walk over a stale view is harmless, being
 * bounded and reading nothing but block ids or NONE, and the CAS that succeeds acquires what the pushers did to
 * the blocks before pushing them.
 * Reached from Router_Init and Router_Push (a sender's refill), and from Router_Init, Router_Progress and
 * Router_Reset whenever the service's stash runs empty.
 */
inline u32 Router_node::free_pop_many(u32 k, u32 *out)
{
	u64 w = free_top.load(std::memory_order_acquire);
	for (;;) {
		u32 n = 0;
		u32 cur = (u32) w;
		while (n < k && cur != ROUTER_NONE) {
			out[n] = cur;
			n += 1;
			cur = (u32) blk_link[cur].load(std::memory_order_relaxed);
		}
		if (n == 0)
			return 0;
		u64 w2 = (((w >> 32) + 1) << 32) | cur;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_acquire))
			return n;
	}
}


/******************************** the sender path ********************************/

/*
 * A sender's cache filled off the free stack, a batch in one CAS, every block zeroed here rather than at its
 * install: lossless waits for at least one block, holding nothing meanwhile, lossy takes what there is.
 * Reached from Router_Init (the sender's constructor) and Router_Push, by the sealer whose install emptied it.
 */
inline void Router_node::refill(Router_thread &s)
{
	u32 n = free_pop_many(ROUTER_BATCH, s.cache);
	while (n == 0 && not lossy) {
		cpu_relax();
		n = free_pop_many(ROUTER_BATCH, s.cache);
	}
	for (u32 i = 0; i < n; i++)
		zero_valid(s.cache[i]);
	s.cache_n = n;
}

/* 
 * Clears a block about to be installed: none of its lines is written(sets the
 * number of lines it)  contains to zero).  Reached from Router_Init and
 * Router_Push (a sender's refill), Router_Progress (a repair, the F-scan's
 * fresh block) and Router_Reset (a repair). 
 */
inline void Router_node::zero_valid(u32 blk)
{
	for (u32 k = 0; k < L; k++)
		n_valid[(size_t) nv_stride * blk + k].store(0, std::memory_order_relaxed);
}

/*
 * Put a sealed block onto the sealed stack, its destination in the link's high word 
 * with the BARE bit when the sealer left the slot without a block
 */
inline void Router_node::seal(int d, u32 blk, bool bare)
{
	u64 hi = ((u64) d << 32) | (bare ? ROUTER_LINK_BARE : 0);
	u32 top = sealed_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[blk].store(hi | top, std::memory_order_relaxed);
		if (sealed_top.compare_exchange_weak(top, blk, std::memory_order_release, std::memory_order_relaxed))
			return;
	}
}

/*
 * Reserve one line slot in the block of destination `d` and write the full line there, so that a block in
 * flight never has a hole; seal the block when the reservation lands on its end, the fresh one coming out of
 * the sender's cache, refilled after the seal once empty; drop the line in lossy mode when the destination has
 * no block or a sealer is installing one: lossy never waits.
 * Reached from Router_Push, by the push that fills the sender's private line for `d`.
 */
inline void Router_node::stage_line(Router_thread &s, int d, const Point *line)
{
	std::atomic<u64> &next = dest[ROUTER_DEST_WORDS * d];
	int n = (int) swc_linesize;
	for (;;) {
		u64 v = next.fetch_add(1ull << 32, std::memory_order_acq_rel);
		u32 k = (u32) (v >> 32);
		u32 blk = (u32) v;
		if (blk == ROUTER_NONE) {         /* the service thread repairs it, once it pops the note */
			s.ctr[ROUTER_DROPPED_SERVICE] += n;
			return;
		}
		if (k < L) {
			Point *pts = (Point *) (pool + (size_t) blk * block_bytes + ROUTER_HDR_BYTES);
			memcpy(pts + (size_t) k * swc_linesize, line, (size_t) n * sizeof(Point));
			n_valid[(size_t) nv_stride * blk + k].store(1, std::memory_order_release);
			return;
		}
		if (k == L) {                     /* the sealer; lossless never finds its cache empty, lossy may */
			u32 fresh = s.cache_n ? s.cache[--s.cache_n] : ROUTER_NONE;
			next.store((u64) fresh, std::memory_order_release);
			seal(d, blk, fresh == ROUTER_NONE);
			if (s.cache_n == 0)
				refill(s);
			continue;
		}
		if (lossy) {                      /* k > L: a sealer is installing, and lossy waits for nobody */
			s.ctr[ROUTER_DROPPED_SERVICE] += n;
			return;
		}
		for (;;) {                        /* at exactly L nobody is installing, and I may become the sealer */
			u64 w = next.load(std::memory_order_acquire);
			if ((u32) (w >> 32) <= L)
				break;
			cpu_relax();
		}
	}
}

/* Sender.  Lossy: wait-free, may drop the line; lossless: may spin, never loses the point. */
inline void Router_Push(u64 a, u64 b, int d, Router_thread &rt)
{
	Point *line = rt.swc + (size_t) d * rt.swc_linesize;
	u64 &count = line[rt.swc_linesize - 1].val;
	u64 n = count;
	line[n].key = a;
	line[n].val = b;                  /* the line's last point lands on the count, read already */
	n += 1;
	rt.ctr[ROUTER_PUSHED] += 1;
	if (n < rt.swc_linesize) {
		count = n;
		return;
	}
	rt.node.stage_line(rt, d, line);
	count = 0;                        /* after the copy: the word held a point until now */
}

/*
 * Sender.  Append every partial line's points to its destination's closing buffer, contiguous with the other
 * senders' thanks to the offset a fetch_add hands out, then announce; nothing more from this sender until
 * Router_Reset.  Nothing partial ever enters a block before the F-scan copies the closing buffers into blocks.
 */
inline void Router_Close(Router_thread &rt)
{
	Router_node &rn = rt.node;
	for (int d = 0; d < rn.F; d++) {
		Point *line = rt.swc + (size_t) d * rt.swc_linesize;
		u64 &count = line[rt.swc_linesize - 1].val;
		if (count == 0)
			continue;
		u64 off = rn.dest[ROUTER_DEST_WORDS * d + 1].fetch_add(count, std::memory_order_relaxed);
		memcpy(rn.partial + (size_t) d * rn.partial_cap + off, line, (size_t) count * sizeof(Point));
		count = 0;
	}
	rt.closed.store(1, std::memory_order_release);
}


/******************************** the receiver's calls ********************************/

/*
 * Receiver.  Up to `buffer_size` points into `buffer` (2 u64 each), out of the blocks the inbox names, in
 * order; returns how many.  A block read through goes onto the free stack, whose CAS is the release that
 * orders these reads before the block's reuse, and is counted after that CAS, which keeps the count behind it.
 */
inline size_t Router_Pop(u64 *buffer, size_t buffer_size, Router_thread &rt)
{
	Router_node &rn = rt.node;
	size_t got = 0;
	while (got < buffer_size) {
		if (rt.cur_blk == ROUTER_NONE) {
			Point e;
			if (not rt.inbox.pop(e))
				break;
			rt.cur_blk = (u32) e.key;
			rt.cur_count = (u32) e.val;
			rt.cur_off = 0;
		}
		size_t take = rt.cur_count - rt.cur_off;
		if (take > buffer_size - got)
			take = buffer_size - got;
		const Point *pts = (const Point *) (rn.pool + (size_t) rt.cur_blk * rn.block_bytes + ROUTER_HDR_BYTES);
		memcpy(buffer + 2 * got, pts + rt.cur_off, take * sizeof(Point));
		got += take;
		rt.cur_off += (u32) take;
		if (rt.cur_off == rt.cur_count) {
			rn.free_push(rt.cur_blk);
			rt.ctr[ROUTER_BLOCKS] += 1;
			rt.cur_blk = ROUTER_NONE;
		}
	}
	rt.ctr[ROUTER_POPPED] += got;
	return got;
}

/* Receiver.  Nothing can arrive any more and nothing is left to read: the flag first, then emptiness, so that
 * a block pushed between the two reads is still popped. */
inline bool Router_Test_drained(Router_thread &rt)
{
	if (rt.node.input_closed.load(std::memory_order_acquire) == 0)
		return false;
	return rt.cur_blk == ROUTER_NONE && rt.inbox.empty();
}

}
#endif
