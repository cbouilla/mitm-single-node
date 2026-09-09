#ifndef MITM_PCS_WALKER
#define MITM_PCS_WALKER

#include <cassert>
#include <algorithm>

#include "tools.hpp"
#include "router/router.hpp"
#include "pcs/params.hpp"
#include "pcs/shared.hpp"
#include "pcs/trail.hpp"

/* PCS's producer: it walks trails, ships their endpoints, and resolves what its dict thread finds. */

namespace mitm::pcs {

/* start chain k from the next chain index j (stepping by jinc) that is not itself a distinguished point */
static inline void start_chain(const Params &params, u64 out_mask, u64 root_seed, u64 &j,
                               u64 x[], u64 len[], u64 seed[], u64 jinc, int k)
{
	u64 start;
	for (;;) {
		j += jinc;
		start = (root_seed + j * params.multiplier) & out_mask;
		if (not is_distinguished_point(start, params.threshold))
			break;
	}
	x[k] = start;
	len[k] = 0;
	seed[k] = j;
}

/*
 * Retire up to `budget` candidates (0 == no limit) from `coll_q`, the queue of this walker's own dict
 * thread.  A step costs a whole vmixf however few slots are in flight, so in the steady state we stop as
 * soon as nothing more is in hand and the batch is not full, and let the rest wait for the next chunk:
 * running the batch down to its last candidate would spend full-width steps on one or two live lanes.
 * `drain` is the end of the round, where every candidate must be retired whatever it costs.  `trail` is
 * the scalar path's recording buffer, and is unused when vlen > 1.
 */
template <class ProblemWrapper>
void service_collisions(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], const Params &params,
                        Shared &shared, CollisionQueue &coll_q, VecResolver<ProblemWrapper> &resolver,
                        size_t budget, bool drain, u64 trail[])
{
	constexpr int vlen = ProblemWrapper::vlen;

	if constexpr (vlen == 1) {
		for (size_t c = 0; resolver.refill(coll_q) > 0; c++) {
			if (budget && c >= budget)
				return;
			CollisionCandidate cand = resolver.pending[resolver.first++];
			resolver.n_pending -= 1;
			assert(cand.i == shared.header.i);
			resolve_collision(wrapper, ctr, hll, params, shared, cand, trail);
		}
		return;
	} else {
		u64 retired = resolver.n_retired;
		for (;;) {
			resolver.fill(wrapper, ctr, params, shared, coll_q);
			if (resolver.n_busy == 0)
				return;                /* nothing in hand, nothing queued, nothing in flight */
			if (not drain && resolver.n_busy < resolver.nslots && resolver.refill(coll_q) == 0)
				return;                /* a partial batch: park it rather than step it at this price */
			resolver.step(wrapper, ctr, hll, params, shared);
			if (budget && resolver.n_retired - retired >= budget)
				return;
		}
	}
}


/*
 * A walker's round: walk vlen trails in lockstep and push every distinguished point to the shard its
 * endpoint names, retiring between chunks the candidates its own dict thread has queued.  Chain indices
 * are strided by the number of walkers from this one's global index, so no two walkers ever share a
 * chain and the partition needs no message.
 *
 * The controller ends the round: `round_over` stops the walking, and what follows is the drain, where
 * every candidate left must be resolved before the next round mixes differently.  A walker may not
 * leave before its dict thread has stopped pushing AND the queue it shares with its fellows is empty.
 */
template <class ProblemWrapper>
void walker_round(Router_thread &rt, const ProblemWrapper &wrapper, const Params &params, Shared &shared,
                  u64 *ctr, u8 *hll, VecResolver<ProblemWrapper> &resolver, u64 *trail, int r)
{
	constexpr int vlen = ProblemWrapper::vlen;
	CollisionQueue &coll_q = *shared.chan[r].q;
	const u64 i = shared.header.i;
	const u64 root_seed = shared.header.root_seed;
	const int jbits = params.jbits;
	const u64 jmask = make_mask(jbits);
	const u64 jinc = Router_num_send(rt);
	const u64 n_recv = Router_num_recv(rt);

	/* state of the vlen chains being walked */
	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 len[vlen];                  /* steps walked so far */
	u64 seed[vlen];                 /* the chain index j of each */

	u64 j = Router_rank(rt);
	for (int k = 0; k < vlen; k++)
		start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, jinc, k);
	assert((j & jmask) == j);

	while (not shared.round_over.load(std::memory_order_relaxed)) {
		service_collisions(wrapper, ctr, hll, params, shared, coll_q, resolver, params.coll_per_chunk,
		                   false, trail);

		for (size_t c = 0; c < params.chunk_size; c++) {
			wrapper.vmixf(i, x, y);

			for (int k = 0; k < vlen; k++) {
				len[k] += 1;
				x[k] = y[k];
				bool dp = is_distinguished_point(x[k], params.threshold);
				bool failure = (len[k] == params.dp_max_it);

				if (dp) {
					ctr[N_DP] += 1;
					ctr[N_POINTS_TRAILS] += len[k];
					u64 l = std::min(len[k], params.len_sat);
					Router_Push(x[k], (seed[k] & jmask) | (l << jbits),
					            (int) (x[k] % n_recv), rt);
				}
				if (dp || failure) {
					if (failure && not dp)
						ctr[BAD_DP] += 1;
					start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, jinc, k);
					assert((j & jmask) == j);
				}
			}
		}
		/* the chunk always runs to completion: exactly one vmixf per turn */
		ctr[N_EVAL] += params.chunk_size * vlen;
	}
	Router_Close(rt);

	/* the drain: `done` first, then the queue, or a run pushed between the two would be left behind */
	for (;;) {
		service_collisions(wrapper, ctr, hll, params, shared, coll_q, resolver, 0, true, trail);
		if (shared.chan[r].done.load(std::memory_order_acquire) && coll_q.is_empty())
			break;
		cpu_relax();
	}
	assert(resolver.n_busy == 0);
	assert(resolver.n_pending == 0);
}

}
#endif
