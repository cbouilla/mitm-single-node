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
 * A located collision (x0, x1), inputs of the same evaluation: tally it, fold it into the node's
 * HyperLogLog, test the pair (PROTOCOL.md §3.4).  Shared by the scalar and the vectorized resolver.
 * True == it was the golden pair, already published.
 */
template<class ProblemWrapper>
bool retire_collision(const ProblemWrapper &wrapper, u64 ctr[], SharedContext &shared,
                      u64 seed0, u64 seed1, u64 x0, u64 x1, u64 len0, u64 len1)
{
	const u64 i = shared.i;
	if (x0 == x1) {
		ctr[BAD_COLLISION] += 1;
		return false;
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

	if (not wrapper.mix_good_pair(i, x0, x1))
		return false;
	printf("\nFound golden collision! i=%" PRIx64 " root_seed=%" PRIx64 " seed0=%" PRIx64
	       ". Dict --> seed1=%" PRIx64 "\n", i, shared.root_seed, seed0, seed1);
	shared.set_golden(i, x0, x1);
	return true;
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
 * right one.  dp_max_it words of stack, and the re-walk stops there rather than overrun them.
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
		if (len1 + 1 >= maxit) {          /* longer than any trail of the round: not the stored one */
			ctr[N_EVAL] += len1;
			ctr[BAD_DP] += 1;
			return nullopt;
		}
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
 * Resolve one candidate, scalar: walk both trails, locate the collision, retire it.  The path taken
 * when the problem has no vector implementation (vlen == 1); VecResolver is the other one.
 */
template<class ProblemWrapper>
void resolve_collision(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
                       SharedContext &shared, u64 seed0, u64 end, u64 len0, u64 seed1, u64 len1_maybe)
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
		return;                 /* robin-hood, or dict false positive */

	auto [x0, x1, len1] = *collision;
	assert(len1_maybe == 0 || len1_maybe == len1);
	retire_collision(wrapper, ctr, shared, seed0, seed1, x0, x1, len0, len1);
}


/******************************* the vectorized resolver **********************/

/*
 * Why: a round spends ~2/beta of its evaluations locating collisions, and a scalar resolver makes
 * each of them vlen times dearer than the ones that produce DPs -- which is most of a walker's time
 * once the dictionary fills.  This one keeps vlen/2 candidates in flight and steps them all with one
 * vmixf.  A candidate walks its two chains through three phases, one lane per chain that is moving:
 *
 *   MEASURE  the stored length had saturated: re-walk chain 1 to its DP to learn it (chain 0 parked)
 *   ALIGN    step the longer chain until both are the same distance from their shared endpoint
 *   MARCH    step both and compare, until they meet or the shorter trail runs out
 *
 * MEASURE and ALIGN hold one lane, MARCH two, so vlen/2 candidates never ask for more than vlen lanes
 * and the pool cannot run dry.  Unlike walk_nolen1() this re-walks chain 1 instead of recording it,
 * which trades a few evaluations for a bounded per-lane footprint.  PROTOCOL.md §3.2.
 */
template<class ProblemWrapper>
struct alignas(sizeof(u64) * ProblemWrapper::vlen) VecResolver {
	static constexpr int vlen = ProblemWrapper::vlen;
	/* vlen == 1 has no vfg(): the scalar resolver runs instead and this object is never touched */
	static constexpr int nslots = (vlen > 1) ? vlen / 2 : 1;
	enum slot_phase {SLOT_EMPTY, SLOT_MEASURE, SLOT_ALIGN, SLOT_MARCH};

	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));   /* the point each lane is walking */
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));   /* one vmixf of x */

	int phase[nslots];        /* what each slot is doing */
	int lane0[nslots];        /* the lane walking chain 0; -1 == parked, still at start0 */
	int lane1[nslots];        /* the lane walking chain 1; -1 == parked, still at start1 */
	u64 seed0[nslots];        /* chain 0: its index, ... */
	u64 start0[nslots];       /* ... where it starts, ... */
	u64 len0[nslots];         /* ... and how long its trail is */
	u64 seed1[nslots];        /* chain 1: its index, ... */
	u64 start1[nslots];       /* ... where it starts, ... */
	u64 len1[nslots];         /* ... and how long its trail is; counted up during MEASURE */
	u64 end0[nslots];         /* the shared endpoint as a dictionary key: MEASURE checks it */
	u64 remaining[nslots];    /* steps left in MEASURE, ALIGN, then MARCH */

	int free_lane[vlen];      /* the lanes no slot is using */
	int n_free;               /* how many of them */
	int n_busy;               /* slots that are not SLOT_EMPTY: the walker is idle when this is 0 */
	u64 n_retired;            /* candidates finished or abandoned, ever: what coll_per_chunk budgets */

	VecResolver()
	{
		for (int l = 0; l < vlen; l++) {
			x[l] = 0;                    /* a free lane still goes through vmixf: keep it in range */
			free_lane[l] = l;
		}
		n_free = vlen;
		n_busy = 0;
		n_retired = 0;
		for (int s = 0; s < nslots; s++) {
			phase[s] = SLOT_EMPTY;
			lane0[s] = -1;
			lane1[s] = -1;
		}
	}

	/* hand a slot's lanes back and free it.  Every path that abandons or completes a candidate ends here */
	void release(int s)
	{
		if (lane0[s] >= 0) {
			x[lane0[s]] = 0;
			free_lane[n_free++] = lane0[s];
			lane0[s] = -1;
		}
		if (lane1[s] >= 0) {
			x[lane1[s]] = 0;
			free_lane[n_free++] = lane1[s];
			lane1[s] = -1;
		}
		phase[s] = SLOT_EMPTY;
		n_busy -= 1;
		n_retired += 1;
	}

	/*
	 * Both chains are now the same distance from their endpoint: give each a lane and step them
	 * together.  Equal points here mean one trail is a suffix of the other, which is no collision.
	 */
	void begin_march(int s, u64 ctr[])
	{
		if (lane0[s] < 0) {
			lane0[s] = free_lane[--n_free];
			x[lane0[s]] = start0[s];
		}
		if (lane1[s] < 0) {
			lane1[s] = free_lane[--n_free];
			x[lane1[s]] = start1[s];
		}
		if (x[lane0[s]] == x[lane1[s]]) {
			ctr[BAD_WALK_ROBINHOOD] += 1;
			release(s);
			return;
		}
		remaining[s] = std::min(len0[s], len1[s]);
		phase[s] = SLOT_MARCH;
	}

	/*
	 * Both trail lengths are known: step the longer chain down to the shorter one's length.  The chain
	 * that does not move stays parked at its start and costs no lane until the march.
	 */
	void begin_align(int s, int held, u64 ctr[])
	{
		assert(held >= 0);                                 /* fill() and MEASURE each hold exactly one */
		lane0[s] = -1;
		lane1[s] = -1;
		if (len0[s] == len1[s]) {
			lane0[s] = held;
			x[held] = start0[s];
			begin_march(s, ctr);
			return;
		}
		if (len0[s] > len1[s]) {
			lane0[s] = held;
			x[held] = start0[s];
			remaining[s] = len0[s] - len1[s];
		} else {
			lane1[s] = held;
			x[held] = start1[s];
			remaining[s] = len1[s] - len0[s];
		}
		phase[s] = SLOT_ALIGN;
	}

	/* start queued candidates in the empty slots.  Returns how many, so a caller can tell an empty queue */
	int fill(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params, SharedContext &shared)
	{
		int started = 0;
		for (int s = 0; s < nslots; s++) {
			if (phase[s] != SLOT_EMPTY)
				continue;
			if (n_free == 0)
				break;                     /* every lane is walking: the rest waits its turn */
			CollisionCandidate c;
			if (not shared.coll_q.pop(c))
				break;                     /* nothing queued */
			assert(c.i == shared.i);
			seed0[s] = c.seed0;
			seed1[s] = c.seed1;
			len0[s] = c.len0;
			len1[s] = c.len1_maybe;
			end0[s] = c.end;
			start0[s] = (shared.root_seed + params.multiplier * c.seed0) & wrapper.out_mask;
			start1[s] = (shared.root_seed + params.multiplier * c.seed1) & wrapper.out_mask;
			int l = free_lane[--n_free];
			n_busy += 1;
			started += 1;
			if (c.len1_maybe == 0) {           /* the stored length saturated: re-walk chain 1 for it */
				lane0[s] = -1;
				lane1[s] = l;
				x[l] = start1[s];
				len1[s] = 0;
				remaining[s] = params.dp_max_it;
				phase[s] = SLOT_MEASURE;
			} else {
				lane0[s] = l;
				lane1[s] = -1;
				x[l] = start0[s];
				phase[s] = SLOT_ALIGN;
				begin_align(s, l, ctr);    /* it picks the chain that has to move, or marches at once */
			}
		}
		return started;
	}

	/*
	 * One vmixf over every lane, then one step of every slot.  N_EVAL counts the evaluations a scalar
	 * resolver would have made, not the lanes burnt, so that it stays comparable across the two paths.
	 */
	void step(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params, SharedContext &shared)
	{
		wrapper.vmixf(shared.i, x, y);

		for (int s = 0; s < nslots; s++) {
			if (phase[s] == SLOT_MEASURE) {
				int l = lane1[s];
				x[l] = y[l];
				len1[s] += 1;
				remaining[s] -= 1;
				ctr[N_EVAL] += 1;
				if (is_distinguished_point(x[l], params.threshold)) {
					if (x[l] / params.n_inserters == end0[s]) {
						begin_align(s, l, ctr);
					} else {
						ctr[BAD_WALK_NONCOLLIDING] += 1;   /* not the trail the dictionary meant */
						release(s);
					}
				} else if (remaining[s] == 0) {
					ctr[BAD_DP] += 1;                  /* the re-walk never reached a DP */
					release(s);
				}
			} else if (phase[s] == SLOT_ALIGN) {
				int l = (lane0[s] >= 0) ? lane0[s] : lane1[s];
				x[l] = y[l];
				remaining[s] -= 1;
				ctr[N_EVAL] += 1;
				if (remaining[s] == 0)
					begin_march(s, ctr);
			} else if (phase[s] == SLOT_MARCH) {
				int a = lane0[s];
				int b = lane1[s];
				ctr[N_EVAL] += 2;
				if (y[a] == y[b]) {
					retire_collision(wrapper, ctr, shared, seed0[s], seed1[s],
					                 x[a], x[b], len0[s], len1[s]);
					release(s);
					continue;
				}
				x[a] = y[a];
				x[b] = y[b];
				remaining[s] -= 1;
				if (remaining[s] == 0) {
					ctr[BAD_WALK_NONCOLLIDING] += 1;   /* the trails never met: a false positive */
					release(s);
				}
			}
		}
	}
};


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


/*
 * Retire up to `budget` collision candidates (0 == no limit).  A step costs a whole vmixf however few
 * slots are in flight, so in the steady state we stop as soon as the queue is empty and the batch is
 * not full, and let the rest wait for the next chunk: running the batch down to its last candidate
 * would spend full-width steps on one or two live lanes.  `drain` is the end of the round, where
 * every candidate must be retired whatever it costs.  PROTOCOL.md §3.2, §3.3.
 */
template <class ProblemWrapper>
void service_collisions(const ProblemWrapper &wrapper, u64 ctr[], const Parameters &params,
                        SharedContext &shared, VecResolver<ProblemWrapper> &resolver,
                        size_t budget, bool drain)
{
	constexpr int vlen = ProblemWrapper::vlen;

	if constexpr (vlen == 1) {
		for (size_t c = 0; not shared.coll_q.is_empty(); c++) {
			if (budget && c >= budget)
				return;
			CollisionCandidate cand;
			if (not shared.coll_q.pop(cand))
				return;
			assert(cand.i == shared.i);
			resolve_collision(wrapper, ctr, params, shared, cand.seed0, cand.end, cand.len0,
			                  cand.seed1, cand.len1_maybe);
		}
		return;
	} else {
		u64 retired = resolver.n_retired;
		for (;;) {
			resolver.fill(wrapper, ctr, params, shared);
			if (resolver.n_busy == 0)
				return;                /* queue empty, nothing in flight */
			if (not drain && resolver.n_busy < resolver.nslots && shared.coll_q.is_empty())
				return;                /* a partial batch: park it rather than step it at this price */
			resolver.step(wrapper, ctr, params, shared);
			if (budget && resolver.n_retired - retired >= budget)
				return;
		}
	}
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

	/* the candidates being resolved, vlen/2 of them at once (PROTOCOL.md §3.2) */
	VecResolver<ProblemWrapper> resolver;

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
			service_collisions(wrapper, ctr, params, shared, resolver, 0, true);
			assert(resolver.n_busy == 0);
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}

		service_collisions(wrapper, ctr, params, shared, resolver, params.coll_per_chunk, false);

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
