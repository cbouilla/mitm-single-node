#ifndef MITM_WALKER
#define MITM_WALKER

#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/*
 * Walking trails: iterate the mixing function to a distinguished point, and turn a
 * dictionary hit into the collision that caused it.
 */

inline bool is_distinguished_point(u64 x, u64 threshold)
{
	return x <= threshold;
}

/*
 * Given two inputs that (maybe) lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * This function assumes that the provided lengths (= distance to the distinguished point)
 * are correct, but it does not assume that the two trails end at the same DP.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
	u64 i, u64 x0, u64 len0, u64 x1, u64 len1__)
{
	/****************************************************************************+
	 *            walk the longest sequence until they are equal                 |
	 * Two chains that leads to the same distinguished point but not necessarily |
	 * have the same length. e.g.                                                |
	 *                                                                           |
	 * chain1: ----------------x-------o                                         |
	 *                        /                                                  |
	 *          chain2: ------                                                   |
	 *                                                                           |
	 * o: is a distinguished point                                               |
	 * x: the collision we're looking for                                        |
	 ****************************************************************************/
	assert(not is_distinguished_point(x1, params.threshold));
	assert(not is_distinguished_point(x0, params.threshold));

	/* move the longest sequence until the remaining number of steps is equal */
	/* to the shortest sequence. */
	u64 len1 = len1__;
	ctr[N_EVAL] += (len0 > len1) ? len0 - len1 : len1 - len0;   /* the two loops below */
	for (; len0 > len1; len0--)
		x0 = wrapper.mixf(i, x0);
	for (; len0 < len1; len1--)
		x1 = wrapper.mixf(i, x1);

	if (x0 == x1) { /* robin-hood */
		ctr[BAD_WALK_ROBINHOOD] += 1;
		return nullopt;
	}

	/* now both sequences needs exactly `len` steps to reach the common distinguished point */
	for (u64 j = 0; j < len0; ++j) {
		/* walk them together and check each time if their output are equal     */
		/* return as soon equality is found. */
		u64 y0 = wrapper.mixf(i, x0);
		u64 y1 = wrapper.mixf(i, x1);
		ctr[N_EVAL] += 2;

		/* First, do the outputs collide? If yes, return true and exit. */
		if (y0 == y1) {
			/* careful: x0 & x1 contain inputs before mixing */
			return optional(tuple(x0, x1, len1__));
		}
		x0 = y0;
		x1 = y1;
	}

	if (x0 != x1)    /* false positive from the dictionnary */
		ctr[BAD_WALK_NONCOLLIDING] += 1;
	return nullopt;
}

/*
 * Given two inputs that (maybe) lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * This function assumes that len0 (= distance to the distinguished point)
 * is correct, but it does not assume that the two trails end at the same DP.
 * `end` is [end of trail] / params.n_inserters
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk_nolen1(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
	u64 i, u64 x0, u64 len0, u64 end0, u64 x1)
{
	/****************************************************************************+
	 *            walk the longest sequence until they are equal                 |
	 * Two chains that leads to the same distinguished point but not necessarily |
	 * have the same length. e.g.                                                |
	 *                                                                           |
	 * chain1: ----------------x-------o                                         |
	 *                        /                                                  |
	 *          chain2: ------                                                   |
	 *                                                                           |
	 * o: is a distinguished point                                               |
	 * x: the collision we're looking for                                        |
	 ****************************************************************************/

	/* the distance from x1 to a distinguished point is unknown.
	 * We need to walk the trail again, but we save all intermediate points.
	 */
	u64 maxit = params.dp_max_it;
	u64 trail1[maxit];
	trail1[0] = x1;
	u64 len1 = 0;
	assert(not is_distinguished_point(x1, params.threshold));
	assert(not is_distinguished_point(x0, params.threshold));
	for (;;) {
		len1 += 1;
		x1 = wrapper.mixf(i, x1);
		trail1[len1] = x1;
		if (is_distinguished_point(x1, params.threshold))
			break;
	}
	ctr[N_EVAL] += len1;                                     /* one per turn of the loop */

	if (x1 / params.n_inserters != end0) {
		ctr[BAD_WALK_NONCOLLIDING] += 1;
		return nullopt;
	}

	/* move the longest sequence until the remaining number of steps is equal */
	/* to the shortest sequence. */
	if (len0 > len1)
		ctr[N_EVAL] += len0 - len1;                          /* the loop below */
	for (; len0 > len1; len0--)
		x0 = wrapper.mixf(i, x0);

	/* at this stage, len0 <= len1 */
	x1 = trail1[len1 - len0];
	if (x0 == x1) { /* robin-hood */
		ctr[BAD_WALK_ROBINHOOD] += 1;
		return nullopt;
	}

	/* now both sequences needs exactly `len0` steps to reach the common distinguished point */
	for (u64 j = len1 - len0;; j++) {
		/* walk them together */
		u64 y0 = wrapper.mixf(i, x0);
		u64 y1 = trail1[j+1];
		ctr[N_EVAL] += 1;
		/* do the outputs collide? If yes, return true and exit. */
		if (y0 == y1) {
			/* careful: x0 & x1 contain inputs before mixing */
			return optional(tuple(x0, x1, len1));
		}
		x0 = y0;
		x1 = y1;
	}
}

/*
 * Walker side of the engine: given a dictionary hit, walk both trails to locate the
 * collision and test whether it is the golden pair.  Returns (i, x0, x1) if it is.
 * Tallies into the walker's `ctr` and into the node's HyperLogLog.
 *
 * This is the expensive half of processing a distinguished point, which is exactly
 * why it does not run on the inserter thread that found the hit.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> resolve_collision(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
                                               SharedContext &shared, u64 seed0, u64 end, u64 len0,
                                               u64 seed1, u64 len1_maybe)
{
	const u64 i = shared.i;
	const u64 root_seed = shared.root_seed;
	u64 start0 = (root_seed + params.multiplier * seed0) & wrapper.out_mask;
	u64 start1 = (root_seed + params.multiplier * seed1) & wrapper.out_mask;

	optional<tuple<u64,u64,u64>> collision;
	if (len1_maybe == 0)
		collision = walk_nolen1(wrapper, ctr, params, i, start0, len0, end, start1);
	else
		collision = walk(wrapper, ctr, params, i, start0, len0, start1, len1_maybe);

	if (not collision)
		return nullopt;         /* robin-hood, or dict false positive */

	auto [x0, x1, len1] = *collision;
	assert(len1_maybe == 0 || len1_maybe == len1);
	if (x0 == x1) {
		ctr[BAD_COLLISION] += 1;
		return nullopt;    /* duh */
	}

	u64 y0 = wrapper.mix(i, x0);
	u64 y1 = wrapper.mix(i, x1);
	assert(wrapper.mixf(i, x0) == wrapper.mixf(i, x1));

	/* the collision's tallies: one more, and its two trail lengths by rank */
	ctr[N_COLLISIONS] += 1;
	if (len0 < len1) {
		ctr[COLLIDING_LEN_MIN] += len0;
		ctr[COLLIDING_LEN_MAX] += len1;
	} else {
		ctr[COLLIDING_LEN_MIN] += len1;
		ctr[COLLIDING_LEN_MAX] += len0;
	}
	shared.found_collision(std::min(y0, y1), std::max(y0, y1));    /* the node's HyperLogLog */

	if (wrapper.mix_good_pair(i, x0, x1)) {
		printf("\nFound golden collision! i=%" PRIx64 " root_seed=%" PRIx64 " seed0=%" PRIx64 ". Dict --> seed1=%" PRIx64 "\n",
			i, root_seed, seed0, seed1);
		return optional(tuple(i, x0, x1));
	}
	return nullopt;
}


static void start_chain(const Parameters &params, u64 out_mask, u64 root_seed, u64 &j,
                        u64 x[], u64 len[], u64 seed[], u64 jinc, int k)
{
	u64 start;
	for (;;) {
		j += jinc;
		start = (root_seed + j * params.multiplier) & out_mask;
		if (not is_distinguished_point(start, params.threshold))  // refuse to start from a DP
			break;
	}
	x[k] = start;
	len[k] = 0;
	seed[k] = j;
}


/*
 * Resolve one queued collision candidate.  Returns false if the queue was empty.
 */
template <class ProblemWrapper>
bool service_collision(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
                       SharedContext &shared)
{
	CollisionCandidate c;
	if (not shared.coll_q.pop(c))
		return false;
	assert(c.i == shared.i);
	auto sol = resolve_collision(wrapper, ctr, params, shared,
	                             c.seed0, c.end, c.len0, c.seed1, c.len1_maybe);
	if (sol) {
		auto [gi, x0, x1] = *sol;
		shared.set_golden(gi, x0, x1);
	}
	return true;
}


/*
 * A walker thread walks `vlen` trails in lockstep and ships every distinguished
 * point it reaches to the comm thread over its own SPSC queue.  Between chunks it
 * also picks up collision candidates that inserter threads have queued, and pays for
 * the expensive part -- walking both trails and testing the pair.
 *
 * The chains are private to this thread, so they are locals; its outgoing queue and
 * its tallies are in its ThreadContext, the collision queue and the HyperLogLog are
 * the node's, in the SharedContext; all are taken by reference up front to keep the
 * hot loop free of the indirection.  The wrapper is stateless and shared read-only
 * by every walker.  `walker_index` is 0-based within
 * this rank -- combined with the rank it gives the global walker index, which seeds
 * the chain counter exactly as `local_rank` did.
 *
 * Winding down is driven by ctx.state (see thread_state in comm.hpp).
 */
template <class ProblemWrapper>
void walker_thread(ThreadContext &ctx, const ProblemWrapper &wrapper, const Parameters &params,
                   SharedContext &shared, int walker_index)
{
	constexpr int vlen = ProblemWrapper::vlen;
	SPSCQueue &out = *ctx.q;
	u64 *ctr = ctx.ctr;
	CollisionQueue &coll_q = shared.coll_q;

	const u64 i = shared.i;
	const u64 root_seed = shared.root_seed;

	int jbits = params.jbits;
	u64 jmask = make_mask(jbits);
	(void) jmask;

	/* state of the vlen chains being walked */
	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 len[vlen], seed[vlen];

	/* same striding scheme as before, with the global walker index in place of local_rank */
	u64 j = (u64) params.rank * params.walkers_per_node + walker_index;
	for (int k = 0; k < vlen; k++)
		start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, params.n_walkers, k);
	assert((j & jmask) == j);

	for (;;) {
		int st = ctx.state.load(std::memory_order_acquire);

		if (st == HOLD) {
			/* the comm thread wants us to stop producing.  Acknowledge, so that it
			   can tell an empty queue from a merely idle-looking one, then carry on
			   with collisions. */
			ctx.state.store(HELD, std::memory_order_release);
			continue;
		}

		if (st == DRAIN) {
			/* every inserter has gone quiet, so no new candidate can appear: finish
			   the queue and this round is over for us */
			while (not coll_q.is_empty())
				service_collision(wrapper, ctr, params, shared);
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}

		/* Retire queued collisions before walking the next chunk. */
		for (size_t c = 0; not coll_q.is_empty(); c++) {
			if (params.coll_per_chunk && c >= params.coll_per_chunk)
				break;
			service_collision(wrapper, ctr, params, shared);
		}

		if (st == HELD) {
			cpu_relax();                     /* held: collisions only, no new points */
			continue;
		}

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
					DP p = {seed[k], x[k], len[k]};
					if (not out.push(p))
						ctr[DROP_WALKERQ] += 1;
				}
				if (dp || failure) {
					if (failure && not dp)
						ctr[BAD_DP] += 1;
					start_chain(params, wrapper.out_mask, root_seed, j,
					            x, len, seed, params.n_walkers, k);
					assert((j & jmask) == j);
				}
			}
		}
		/* the chunk always runs to completion: exactly one vmixf per turn */
		ctr[N_EVAL] += params.chunk_size * vlen;
	}
}

}
#endif
