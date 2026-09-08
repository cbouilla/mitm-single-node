#ifndef MITM_ROUTER_WORKERS
#define MITM_ROUTER_WORKERS

#include "common.hpp"

namespace mitm {

/******************************** the free ring ********************************/

/*
 * Put n blocks onto the free ring, a run of tickets then their cells, each once the popper of its previous
 * element has stored it free -- at most that popper's last stores away, the ring being never full.  The cell's
 * store releases everything the pusher did to the block.
 * Reached from Router_Pop (one block) and from spill, in Router_Init and Router_Progress (a batch).
 */
inline void Router_node::free_push(const u32 *blks, u32 n)
{
	u64 t = free_in.fetch_add(n, std::memory_order_relaxed);
	for (u32 i = 0; i < n; i++) {
		std::atomic<u64> &cell = free_cell[(t + i) & free_mask];
		while ((u32) (cell.load(std::memory_order_acquire) >> 32) != (u32) (t + i))
			cpu_relax();
		cell.store(((u64) (u32) (t + i + 1) << 32) | blks[i], std::memory_order_release);
	}
}

/*
 * Pop up to k blocks off the free ring: the ready prefix at free_out, claimed whole by one CAS; zero when the
 * element at free_out is not there (nothing free, or a push under way), never a wait.  A ready cell changes only
 * when claimed, and a claim moves free_out, which fails the CAS: the ids read before it are the blocks.  The
 * acquire loads pair with the pushers' release stores: what they did to a block is visible to whoever pops it.
 * Reached from Router_Push (a sender's refill) and from the service whenever its stash runs empty.
 */
inline u32 Router_node::free_pop_many(u32 k, u32 *out)
{
	u64 pos = free_out.load(std::memory_order_relaxed);
	u32 n = 0;
	while (n == 0) {
		u64 c = free_cell[pos & free_mask].load(std::memory_order_acquire);
		int32_t dif = (int32_t) ((u32) (c >> 32) - (u32) (pos + 1));
		if (dif < 0)
			return 0;                     /* element pos is not there yet */
		if (dif > 0) {                    /* claimed already: my view of free_out is stale */
			pos = free_out.load(std::memory_order_relaxed);
			continue;
		}
		out[0] = (u32) c;
		n = 1;
		while (n < k) {                   /* the ready prefix behind it */
			c = free_cell[(pos + n) & free_mask].load(std::memory_order_acquire);
			if ((u32) (c >> 32) != (u32) (pos + n + 1))
				break;
			out[n] = (u32) c;
			n += 1;
		}
		if (not free_out.compare_exchange_weak(pos, pos + n, std::memory_order_relaxed))
			n = 0;                        /* lost the race; pos now holds the value that won */
	}
	for (u32 i = 0; i < n; i++)           /* the cells, free for their next round */
		free_cell[(pos + i) & free_mask].store(((u64) (u32) (pos + i + free_mask + 1) << 32) | ROUTER_NONE,
		                                      std::memory_order_release);
	return n;
}


/******************************** the sender path ********************************/

/*
 * A sender's cache filled off the free ring, a batch in one CAS, every block zeroed here rather than at its
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
 * order; returns how many.  A block read through goes onto the free ring, whose cell store is the release that
 * orders these reads before the block's reuse, and is counted after that store, which keeps the count behind it.
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
			rn.free_push(&rt.cur_blk, 1);
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
