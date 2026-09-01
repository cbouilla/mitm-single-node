#ifndef MITM_WALKER
#define MITM_WALKER

#include <cmath>
#include <vector>

#include "trail.hpp"
#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

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
 * A named function rather than a lambda inside walker_thread(), so everything it
 * touches shows up in the signature.
 */
template <class ProblemWrapper>
bool service_collision(ThreadContext &ctx, ProblemWrapper &wrapper, const Parameters &params,
                       RoundState &round, CollisionQueue &coll_q, u64 i, u64 root_seed)
{
	CollisionCandidate c;
	if (not coll_q.pop(c))
		return false;
	/* same invariant: the queue starts every round empty, so a candidate can only
	   ever belong to the round we are in */
	assert(c.i == i);
	auto sol = resolve_collision(wrapper, ctx.ctr, params, i, root_seed,
	                             c.seed0, c.end, c.len0, c.seed1, c.len1_maybe);
	if (sol) {
		auto [gi, x0, x1] = *sol;
		round.set_golden(gi, x0, x1);
	}
	return true;
}


/*
 * A walker thread walks `vlen` trails in lockstep and ships every distinguished
 * point it reaches to the comm thread over its own SPSC queue.  Between chunks it
 * also picks up collision candidates that inserter threads have queued, and pays for
 * the expensive part -- walking both trails and testing the pair.
 *
 * The chains, the problem wrapper and the outgoing queue are private to this thread,
 * so they are locals and arguments; only the tallies the comm thread reads live in
 * ctx.  `walker_index` is 0-based within this rank -- combined with the rank it gives
 * the global walker index, which seeds the chain counter exactly as `local_rank` did.
 *
 * Winding down is driven by ctx.state (see thread_state in comm.hpp).
 */
template <class ProblemWrapper>
void walker_thread(ThreadContext &ctx, const ProblemWrapper &master, const Parameters &params,
                   RoundState &round, SPSCQueue &out, CollisionQueue &coll_q, int walker_index)
{
	constexpr int vlen = ProblemWrapper::vlen;

	/* our own copy: wrapper.n_eval is a per-thread tally */
	ProblemWrapper wrapper = master;
	wrapper.n_eval = 0;

	const u64 i = round.i;
	const u64 root_seed = round.root_seed;

	int jbits = std::log2(10 * params.w) + 8;
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

	u64 n_dp_local = 0;

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
			while (coll_q.count.load(std::memory_order_relaxed) > 0)
				service_collision(ctx, wrapper, params, round, coll_q, i, root_seed);
			ctx.n_dp.store(n_dp_local, std::memory_order_relaxed);
			ctx.n_eval.store(wrapper.n_eval, std::memory_order_relaxed);
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}

		/*
		 * Retire queued collisions before walking the next chunk.  Draining rather
		 * than taking a single one is what keeps up: a chunk produces chunk_size*vlen
		 * points, hence far more than one dictionary hit, so retiring one at a time
		 * loses almost every collision to the queue.  Draining is also self-balancing
		 * -- if the inserters are finding hits faster than we retire them, we spend
		 * more time here and produce fewer points.
		 */
		for (size_t c = 0; coll_q.count.load(std::memory_order_relaxed) > 0; c++) {
			if (params.coll_per_chunk && c >= params.coll_per_chunk)
				break;
			service_collision(ctx, wrapper, params, round, coll_q, i, root_seed);
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
					n_dp_local += 1;
					ctx.ctr.n_points_trails += len[k];
					DP p = {seed[k], x[k], len[k]};
					if (not out.push(p))
						ctx.n_drop_walkerq.fetch_add(1, std::memory_order_relaxed);
				}
				if (dp || failure) {
					if (failure && not dp)
						ctx.ctr.bad_dp += 1;
					start_chain(params, wrapper.out_mask, root_seed, j,
					            x, len, seed, params.n_walkers, k);
					assert((j & jmask) == j);
				}
			}
		}

		/* publish progress for the comm thread's periodic report */
		ctx.n_dp.store(n_dp_local, std::memory_order_relaxed);
		ctx.n_eval.store(wrapper.n_eval, std::memory_order_relaxed);
	}
}

}
#endif
