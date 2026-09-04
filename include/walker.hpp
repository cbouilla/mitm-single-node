#ifndef MITM_WALKER
#define MITM_WALKER

#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/* The walker side: trails, the collision behind a dictionary hit, and the walker thread. */

inline bool is_distinguished_point(u64 x, u64 threshold)
{
	return x <= threshold;
}

/*
 * The collision behind a dictionary hit, both trail lengths known: (x0, x1, len1) with
 * mixf(x0) == mixf(x1), or nothing when one trail is a suffix of the other (robin-hood) or they never
 * meet (a false positive).  Trusts the lengths, not that the trails end at the same DP.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
	u64 i, u64 x0, u64 len0, u64 x1, u64 len1__)
{
	assert(not is_distinguished_point(x1, params.threshold));
	assert(not is_distinguished_point(x0, params.threshold));

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

	for (u64 j = 0; j < len0; ++j) {
		u64 y0 = wrapper.mixf(i, x0);
		u64 y1 = wrapper.mixf(i, x1);
		ctr[N_EVAL] += 2;
		if (y0 == y1) {
			return optional(tuple(x0, x1, len1__));   /* x0, x1: the inputs of the colliding evaluation */
		}
		x0 = y0;
		x1 = y1;
	}

	if (x0 != x1)    /* false positive from the dictionary */
		ctr[BAD_WALK_NONCOLLIDING] += 1;
	return nullopt;
}

/*
 * walk() when the second trail's length is unknown (it saturated in the dictionary): that trail is
 * re-walked and kept, and `end0` (its DP as dictionary key, x / n_inserters) says whether it is the
 * right one.  dp_max_it words of stack.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk_nolen1(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
	u64 i, u64 x0, u64 len0, u64 end0, u64 x1)
{
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

	if (len0 > len1)
		ctr[N_EVAL] += len0 - len1;                          /* the loop below */
	for (; len0 > len1; len0--)
		x0 = wrapper.mixf(i, x0);

	x1 = trail1[len1 - len0];
	if (x0 == x1) { /* robin-hood */
		ctr[BAD_WALK_ROBINHOOD] += 1;
		return nullopt;
	}

	for (u64 j = len1 - len0;; j++) {
		u64 y0 = wrapper.mixf(i, x0);
		u64 y1 = trail1[j+1];
		ctr[N_EVAL] += 1;
		if (y0 == y1) {
			return optional(tuple(x0, x1, len1));   /* x0, x1: the inputs of the colliding evaluation */
		}
		x0 = y0;
		x1 = y1;
	}
}

/*
 * Resolve a candidate: walk both trails, locate the collision, tally it (counters and the node's
 * HyperLogLog), test it.  Returns (i, x0, x1) if it is the golden pair.  The expensive half of a DP,
 * which is why it runs on a walker and not on the inserter that found the hit.
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
		return nullopt;
	}

	u64 y0 = wrapper.mix(i, x0);
	u64 y1 = wrapper.mix(i, x1);
	assert(wrapper.mixf(i, x0) == wrapper.mixf(i, x1));

	ctr[N_COLLISIONS] += 1;
	if (len0 < len1) {
		ctr[COLLIDING_LEN_MIN] += len0;
		ctr[COLLIDING_LEN_MAX] += len1;
	} else {
		ctr[COLLIDING_LEN_MIN] += len1;
		ctr[COLLIDING_LEN_MAX] += len0;
	}
	shared.found_collision(std::min(y0, y1), std::max(y0, y1));

	if (wrapper.mix_good_pair(i, x0, x1)) {
		printf("\nFound golden collision! i=%" PRIx64 " root_seed=%" PRIx64 " seed0=%" PRIx64
		       ". Dict --> seed1=%" PRIx64 "\n", i, root_seed, seed0, seed1);
		return optional(tuple(i, x0, x1));
	}
	return nullopt;
}


/* start chain k from the next chain index j (stepping by jinc) that is not itself a DP.  theta == 1: never returns */
static void start_chain(const Parameters &params, u64 out_mask, u64 root_seed, u64 &j,
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


/* pop and resolve one collision candidate.  False if the queue was empty */
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
 * A walker thread: walks vlen trails in lockstep, ships every DP to the comm thread over its SPSC
 * queue, and between chunks retires the candidates the inserters queued.  Chain indices are strided
 * by n_walkers from the global walker index.  Wind-down: ctx.state, PROTOCOL.md §3.3.
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
	u64 len[vlen];                  /* steps walked so far */
	u64 seed[vlen];                 /* the chain index j of each */

	u64 j = (u64) params.rank * params.walkers_per_node + walker_index;
	for (int k = 0; k < vlen; k++)
		start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, params.n_walkers, k);
	assert((j & jmask) == j);

	for (;;) {
		int st = ctx.state.load(std::memory_order_acquire);

		if (st == HOLD) {
			ctx.state.store(HELD, std::memory_order_release);
			continue;
		}

		if (st == DRAIN) {
			while (not coll_q.is_empty())
				service_collision(wrapper, ctr, params, shared);
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}

		for (size_t c = 0; not coll_q.is_empty(); c++) {
			if (params.coll_per_chunk && c >= params.coll_per_chunk)
				break;
			service_collision(wrapper, ctr, params, shared);
		}

		if (st == HELD) {
			cpu_relax();
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
