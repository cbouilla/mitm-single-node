#ifndef MITM_INSERTER
#define MITM_INSERTER

#include <cassert>

#include "tools.hpp"
#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/******************************* dictionary shard *****************************/

/*
 * This dictionary, when probed with the distinguished point at the end of a trail,
 * should provide (if any) the start (and the length) of another distinguished point
 * that has the same end.
 *
 * One instance per inserter thread, in SharedContext::shards, built by that thread
 * once it is pinned (see run() in engine.hpp): the zero-fill in the constructor is
 * the first write to every page of A, and under Linux's first-touch policy that is
 * what puts the shard on the NUMA node of the CPU that will probe it.  comm.hpp only
 * forward-declares the class: the table holds pointers, and the one SharedContext is
 * built and destroyed in run(), where the class is complete.
 */
class PcsDict {
public:
	u64 jbits, lbits;
	u64 jmask, lmask;
	u64 key_mask;
	const u64 n_slots;     /* size of A */

	vector<u64> A;         // A[i][0:jbits] == j.  A[i][jbits:lbits] == len1.  A[lbits:64] == key bits

	PcsDict(u64 jbits, u64 w) : jbits(jbits), n_slots(w)
	{
		assert(jbits <= 56);
		jmask = make_mask(jbits);
		lmask = make_mask(8);
		lbits = jbits + 8;
		key_mask = (lbits == 64) ? 0 : 0xffffffffffffffff << lbits;
		A.resize(n_slots);     /* zero-filled by the calling thread: the first touch of every page */
	}

	/*
	 * Between two rounds, by the owner alone.  Placement is settled by then, so any
	 * thread could help; those on the shard's NUMA node should (the flush is
	 * bandwidth-bound and worth sharing when the inserters are few).  Which threads
	 * share the shard's NUMA node is known (ThreadContext::numa_node); the collective
	 * flush itself is future work.
	 */
	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i] = 0;
	}

	// return (start', len'), maybe. Return len' == 0 if len' is unknown (has been truncated)
	optional<pair<u64, u64>> pop_insert(u64 end, u64 start, u64 len0)
	{
		u64 idx = end % n_slots;
		u64 key = (end / n_slots) << lbits;

		// first save the previous point
		u64 e = A[idx];
		u64 ekey = e & key_mask;
		u64 elen = (e >> jbits) & lmask;

		/*
		 * heuristic modification from the original algorithm:
		 * we only overwrite an existing point if we have a longer tail.
		 * the alternative is to always insert and forget about it.
		 */
		if (e == 0 || len0 >= elen) {
			if (len0 > lmask)      // saturate the length because be don't use many bits for it.
				len0 = lmask;
			A[idx] = start ^ (len0 << jbits) ^ key;
		}

		if (ekey != key || e == 0)
			return nullopt;

		if (elen == lmask)
			elen = 0;

		return optional(pair(e & jmask, elen));
	}
};

/******************************* the inserter thread **************************/

/*
 * REMARK:
 * Each probe is a random read-modify-write into a multi-GB array: an L3
 * and TLB miss, and they run back-to-back with nothing to overlap them.
 * Software-pipelining is the obvious idea -- issue __builtin_prefetch for
 * a new point, then probe one queued a dozen points earlier, so that many
 * misses are outstanding at once.  That was built and measured: about 5%
 * on a laptop with a 256 MB dictionary, which did not pay for the ring
 * buffer, so it was taken back out.  Worth retrying on cluster hardware
 * (bigger caches, more memory channels) before dismissing it.
 */


/*
 * An inserter thread owns one shard of the distributed dictionary and does exactly
 * one thing with each incoming distinguished point: probe it.  Everything expensive
 * that follows a hit (walking the two trails) is handed to a walker thread through
 * the collision queue.
 *
 * The shard is shared.shards[inserter_index], built by this thread once it was
 * pinned and this thread's alone to probe and flush; the incoming queue and the
 * tallies are in its ThreadContext, the collision queue is the node's.  All are
 * taken by reference up front to keep the loop free of the indirection.  Winding
 * down is driven by ctx.state (see thread_state in comm.hpp).
 */
inline void inserter_thread(ThreadContext &ctx, const Parameters &params, SharedContext &shared,
                            int inserter_index)
{
	SPSCQueue &in = *ctx.q;
	PcsDict &dict = *shared.shards[inserter_index];
	CollisionQueue &coll_q = shared.coll_q;
	u64 *ctr = ctx.ctr;
	static const size_t BATCH = 64;
	DP staging[BATCH];

	for (;;) {
		while (not in.empty()) {
			size_t k = in.pop_bulk(staging, BATCH);
			for (size_t t = 0; t < k; t++) {
				const DP &p = staging[t];
				/* by construction this point belongs to THIS thread's shard */
				u64 key = p.x / params.n_inserters;
				ctr[N_PROBE] += 1;
				auto hit = dict.pop_insert(key, p.seed, p.len);
				if (not hit) {
					/* the trails do not coalesce, just a hash collision */
					ctr[BAD_PROBE] += 1;
					continue;
				}
				auto [seed1, len1_maybe] = *hit;
				CollisionCandidate c = {shared.i, p.seed, p.len, key, seed1, len1_maybe};
				if (not coll_q.push(c))
					ctr[DROP_COLL] += 1;
			}
		}

		// FIRST check state, THEN check emptiness (otherwise the queue could refill between both checks)
		if (ctx.state.load(std::memory_order_acquire) == DRAIN) {
			if (in.empty()) {
				ctx.state.store(QUIESCENT, std::memory_order_release);
				return;
			}
			continue;    // queue is not empty, don't relax
		}

		cpu_relax();
	}
}

}
#endif
