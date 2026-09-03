#ifndef MITM_INSERTER
#define MITM_INSERTER

#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

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
 * The shard, the incoming queue and the tallies are this thread's own, built with
 * its ThreadContext once it was pinned; all three are taken by reference up front to
 * keep the loop free of the indirection.  Winding down is driven by ctx.state (see
 * thread_state in comm.hpp).
 */
inline void inserter_thread(ThreadContext &ctx, const Parameters &params, RoundState &round,
                            CollisionQueue &coll_q)
{
	SPSCQueue &in = *ctx.q;
	PcsDict &dict = *ctx.dict;
	Counters &ctr = ctx.ctr;
	static const size_t BATCH = 64;
	DP staging[BATCH];

	for (;;) {
		while (not in.empty()) {
			size_t k = in.pop_bulk(staging, BATCH);
			for (size_t t = 0; t < k; t++) {
				const DP &p = staging[t];
				/* by construction this point belongs to THIS thread's shard */
				u64 key = p.x / params.n_inserters;
				ctr.c[N_PROBE] += 1;
				auto hit = dict.pop_insert(key, p.seed, p.len);
				if (not hit) {
					/* the trails do not coalesce, just a hash collision */
					ctr.c[BAD_PROBE] += 1;
					continue;
				}
				auto [seed1, len1_maybe] = *hit;
				CollisionCandidate c = {round.i, p.seed, p.len, key, seed1, len1_maybe};
				if (not coll_q.push(c))
					ctr.c[DROP_COLL] += 1;
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
