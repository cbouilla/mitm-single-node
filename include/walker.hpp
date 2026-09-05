#ifndef MITM_WALKER
#define MITM_WALKER

#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

#include "parameters.hpp"
#include "pcs_common.hpp"

namespace mitm::pcs {

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
bool retire_collision(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], SharedContext<Scheme> &shared,
                      u64 seed0, u64 seed1, u64 x0, u64 x1, u64 len0, u64 len1)
{
	const u64 i = shared.header.i;
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
	hll_record(hll, std::min(y0, y1), std::max(y0, y1));

	if (not wrapper.mix_good_pair(i, x0, x1))
		return false;
	printf("\nFound golden collision! i=%" PRIx64 " root_seed=%" PRIx64 " seed0=%" PRIx64
	       ". Dict --> seed1=%" PRIx64 "\n", i, shared.header.root_seed, seed0, seed1);
	shared.set_golden(i, x0, x1);
	return true;
}

/*
 * The collision behind a dictionary hit, both trail lengths known: (x0, x1) with mixf(x0) == mixf(x1),
 * or nothing when one trail is a suffix of the other (robin-hood) or they never meet (a false
 * positive).  Trusts the lengths, not that the trails end at the same DP.
 */
template<class ProblemWrapper>
optional<pair<u64,u64>> walk(const ProblemWrapper &wrapper, u64 ctr[], const Params &params,
	u64 i, u64 x0, u64 len0, u64 x1, u64 len1)
{
	assert(not is_distinguished_point(x1, params.threshold));
	assert(not is_distinguished_point(x0, params.threshold));

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
			return optional(pair(x0, x1));   /* x0, x1: the inputs of the colliding evaluation */
		}
		x0 = y0;
		x1 = y1;
	}

	if (x0 != x1)    /* false positive from the dictionary */
		ctr[BAD_WALK_NONCOLLIDING] += 1;
	return nullopt;
}

/*
 * The trail length of one chain of a candidate, recovered by re-walking it to its distinguished point:
 * what a length that saturated on the wire or in the dictionary costs (PROTOCOL.md §3.2).  Nothing when
 * the chain reaches no distinguished point within dp_max_it steps, or reaches one other than `end`.
 */
template<class ProblemWrapper>
optional<u64> measure_trail(const ProblemWrapper &wrapper, u64 ctr[], const Params &params,
	u64 i, u64 x, u64 end)
{
	assert(not is_distinguished_point(x, params.threshold));
	ctr[N_MEASURE] += 1;
	for (u64 len = 1; len <= params.dp_max_it; len++) {
		x = wrapper.mixf(i, x);
		ctr[N_EVAL] += 1;
		if (not is_distinguished_point(x, params.threshold))
			continue;
		if (x != end) {
			ctr[BAD_WALK_NONCOLLIDING] += 1;      /* not the trail the dictionary meant */
			return nullopt;
		}
		return optional(len);
	}
	ctr[BAD_DP] += 1;
	return nullopt;
}

/*
 * walk() when one of the two trail lengths is unknown: that trail is re-walked into `trail` and kept,
 * which recovers its length and makes the march one evaluation per step instead of two.  The other
 * trail is stepped and its length is trusted; either chain of the candidate can play either role.
 * Returns the inputs of the colliding evaluation, the stepped chain's first, then the recorded trail's
 * length.  `trail` holds dp_max_it + 1 points, which is every trail a round can produce.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk_recorded(const ProblemWrapper &wrapper, u64 ctr[],
	const Params &params, u64 i, u64 x_step, u64 len_step, u64 x_rec, u64 end, u64 trail[])
{
	assert(not is_distinguished_point(x_step, params.threshold));
	assert(not is_distinguished_point(x_rec, params.threshold));
	assert(len_step > 0);

	ctr[N_MEASURE] += 1;
	trail[0] = x_rec;
	u64 len_rec = 0;
	for (;;) {
		len_rec += 1;
		x_rec = wrapper.mixf(i, x_rec);
		trail[len_rec] = x_rec;
		if (is_distinguished_point(x_rec, params.threshold))
			break;
		if (len_rec == params.dp_max_it) {        /* longer than any trail of the round: not the stored one */
			ctr[N_EVAL] += len_rec;
			ctr[BAD_DP] += 1;
			return nullopt;
		}
	}
	ctr[N_EVAL] += len_rec;                       /* one per turn of the loop */

	if (x_rec != end) {
		ctr[BAD_WALK_NONCOLLIDING] += 1;
		return nullopt;
	}

	if (len_step > len_rec)
		ctr[N_EVAL] += len_step - len_rec;        /* the loop below */
	for (; len_step > len_rec; len_step--)
		x_step = wrapper.mixf(i, x_step);

	u64 j = len_rec - len_step;
	if (x_step == trail[j]) { /* robin-hood */
		ctr[BAD_WALK_ROBINHOOD] += 1;
		return nullopt;
	}

	for (; j < len_rec; j++) {
		u64 y_step = wrapper.mixf(i, x_step);
		ctr[N_EVAL] += 1;
		if (y_step == trail[j + 1])
			return optional(tuple(x_step, trail[j], len_rec));
		x_step = y_step;
	}

	ctr[BAD_WALK_NONCOLLIDING] += 1;              /* the trails never met: a false positive */
	return nullopt;
}

/*
 * Resolve one candidate, scalar: recover whichever of the two trail lengths saturated, walk both
 * trails, locate the collision, retire it (PROTOCOL.md §3.2).  The path taken when the problem has no
 * vector implementation (vlen == 1); VecResolver is the other one.  A chain whose length is unknown is
 * the one recorded, and when both are unknown the other is measured first, because walk_recorded()
 * steps it and has to know how far.
 */
template<class ProblemWrapper>
void resolve_collision(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], const Params &params,
                       SharedContext<Scheme> &shared, const CollisionCandidate &c, u64 trail[])
{
	const u64 i = shared.header.i;
	const u64 root_seed = shared.header.root_seed;
	u64 start0 = (root_seed + params.multiplier * c.seed0) & wrapper.out_mask;
	u64 start1 = (root_seed + params.multiplier * c.seed1) & wrapper.out_mask;
	u64 len0 = c.len0_maybe;
	u64 len1 = c.len1_maybe;

	if (len0 != 0 && len1 != 0) {
		auto collision = walk(wrapper, ctr, params, i, start0, len0, start1, len1);
		if (not collision)
			return;                 /* robin-hood, or dict false positive */
		auto [x0, x1] = *collision;
		retire_collision(wrapper, ctr, hll, shared, c.seed0, c.seed1, x0, x1, len0, len1);
		return;
	}

	if (len0 == 0 && len1 == 0) {
		auto measured = measure_trail(wrapper, ctr, params, i, start0, c.end);
		if (not measured)
			return;
		len0 = *measured;
		assert(len0 >= params.len_sat);            /* the wire said "at least that long" */
	}

	if (len1 == 0) {
		auto collision = walk_recorded(wrapper, ctr, params, i, start0, len0, start1, c.end, trail);
		if (not collision)
			return;
		auto [x0, x1, len_rec] = *collision;
		retire_collision(wrapper, ctr, hll, shared, c.seed0, c.seed1, x0, x1, len0, len_rec);
	} else {
		auto collision = walk_recorded(wrapper, ctr, params, i, start1, len1, start0, c.end, trail);
		if (not collision)
			return;
		auto [x1, x0, len_rec] = *collision;
		assert(len_rec >= params.len_sat);         /* likewise: chain 0's length came off the wire */
		retire_collision(wrapper, ctr, hll, shared, c.seed0, c.seed1, x0, x1, len_rec, len1);
	}
}

/******************************* the vectorized resolver **********************/

/*
 * Why: a round spends ~2/beta of its evaluations locating collisions, and a scalar resolver makes
 * each of them vlen times dearer than the ones that produce DPs -- which is most of a walker's time
 * once the dictionary fills.  This one keeps vlen/2 candidates in flight and steps them all with one
 * vmixf.  A candidate walks its two chains through three phases:
 *
 *   MEASURE  a length saturated: re-walk that chain to its DP to learn it -- either chain, or both
 *   ALIGN    step the longer chain until both are the same distance from their shared endpoint
 *   MARCH    step both and compare, until they meet or the shorter trail runs out
 *
 * A busy slot owns two lanes from fill() to release(), so vlen/2 candidates never ask for more than
 * vlen lanes and the pool cannot run dry: an empty slot always finds its pair.  A chain that is not
 * moving is still stepped by vmixf, its lane simply not committed.  Unlike walk_recorded() this
 * re-walks a chain instead of recording it, which trades a few evaluations for a bounded per-lane
 * footprint.  PROTOCOL.md §3.2.
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
	int lane0[nslots];        /* the lane chain 0 owns while the slot is busy; -1 == the slot is empty */
	int lane1[nslots];        /* likewise chain 1 */
	u64 seed0[nslots];        /* chain 0: its index, ... */
	u64 start0[nslots];       /* ... where it starts, ... */
	u64 len0[nslots];         /* ... and how long its trail is; counted up while it measures */
	u64 seed1[nslots];        /* chain 1: its index, ... */
	u64 start1[nslots];       /* ... where it starts, ... */
	u64 len1[nslots];         /* ... and how long its trail is; counted up while it measures */
	bool measuring0[nslots];  /* chain 0 is still walking to its DP: its length saturated on the wire */
	bool measuring1[nslots];  /* chain 1 is still walking to its DP: it saturated in the dictionary */
	u64 end[nslots];          /* the shared endpoint, in full: what a measured chain must reach */
	u64 remaining[nslots];    /* steps left in ALIGN, then in MARCH */

	int free_lane[vlen];      /* the lanes no slot is using */
	int n_free;               /* how many of them */
	int n_busy;               /* slots that are not SLOT_EMPTY: the walker is idle when this is 0 */
	u64 n_retired;            /* candidates finished or abandoned, ever: what coll_per_chunk budgets */

	/* the run this walker took from its group's queue; the scalar resolver draws from it too (§3.2) */
	static constexpr int PENDING = 64;
	CollisionCandidate pending[PENDING];
	int first;                /* the next candidate to start, pending[first] */
	int n_pending;            /* how many are left in the run */

	VecResolver()
	{
		first = 0;
		n_pending = 0;
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
			measuring0[s] = false;
			measuring1[s] = false;
		}
	}

	/* how many candidates are in hand, taking a fresh run from the queue when the last one ran out */
	int refill(CollisionQueue &coll_q)
	{
		if (n_pending > 0)
			return n_pending;
		if (coll_q.is_empty())
			return 0;              /* relaxed read: an idle queue costs no lock traffic (PROTOCOL.md §3.2) */
		first = 0;
		n_pending = (int) coll_q.pop_bulk(pending, PENDING);
		return n_pending;
	}

	/* hand a slot's two lanes back and free it.  Every path that abandons or completes a candidate ends here */
	void release(int s)
	{
		assert(lane0[s] >= 0 && lane1[s] >= 0);
		x[lane0[s]] = 0;
		free_lane[n_free++] = lane0[s];
		lane0[s] = -1;
		x[lane1[s]] = 0;
		free_lane[n_free++] = lane1[s];
		lane1[s] = -1;
		measuring0[s] = false;
		measuring1[s] = false;
		phase[s] = SLOT_EMPTY;
		n_busy -= 1;
		n_retired += 1;
	}

	/*
	 * Both chains are now the same distance from their endpoint: step them together.  Equal points here
	 * mean one trail is a suffix of the other, which is no collision.
	 */
	void begin_march(int s, u64 ctr[])
	{
		if (x[lane0[s]] == x[lane1[s]]) {
			ctr[BAD_WALK_ROBINHOOD] += 1;
			release(s);
			return;
		}
		remaining[s] = std::min(len0[s], len1[s]);
		phase[s] = SLOT_MARCH;
	}

	/*
	 * Both trail lengths are known: rewind both chains and step the longer one down to the shorter
	 * one's length.  The chain that does not move waits at its start, in its own lane.
	 */
	void begin_align(int s, u64 ctr[])
	{
		assert(lane0[s] >= 0 && lane1[s] >= 0);   /* a busy slot owns both lanes, fill() to release() */
		assert(len0[s] > 0 && len1[s] > 0);       /* every length is exact by now: MEASURE is over */
		x[lane0[s]] = start0[s];                  /* a measured chain sits on its endpoint: rewind it */
		x[lane1[s]] = start1[s];
		if (len0[s] == len1[s]) {
			begin_march(s, ctr);
			return;
		}
		remaining[s] = (len0[s] > len1[s]) ? len0[s] - len1[s] : len1[s] - len0[s];
		phase[s] = SLOT_ALIGN;
	}

	/*
	 * One step of a chain re-walking to its distinguished point to recover its length.  False == the
	 * candidate is gone and the slot has been released, so the caller must not touch it again.
	 */
	bool measure_chain(int s, int l, u64 &len, bool &measuring, const Params &params, u64 ctr[])
	{
		x[l] = y[l];
		len += 1;
		ctr[N_EVAL] += 1;
		if (is_distinguished_point(x[l], params.threshold)) {
			if (x[l] != end[s]) {
				ctr[BAD_WALK_NONCOLLIDING] += 1;   /* not the trail the dictionary meant */
				release(s);
				return false;
			}
			measuring = false;
			return true;
		}
		if (len == params.dp_max_it) {
			ctr[BAD_DP] += 1;                  /* the re-walk never reached a DP */
			release(s);
			return false;
		}
		return true;
	}

	/* start queued candidates in the empty slots.  Returns how many, so a caller can tell an empty queue */
	int fill(const ProblemWrapper &wrapper, u64 ctr[], const Params &params, SharedContext<Scheme> &shared,
	         CollisionQueue &coll_q)
	{
		/* the lane budget: two lanes per busy slot, so an empty slot always finds its pair */
		static_assert(vlen == 1 || 2 * nslots <= vlen, "a busy slot owns two of the vlen lanes");
		assert(n_free == vlen - 2 * n_busy);
		int started = 0;
		for (int s = 0; s < nslots; s++) {
			if (phase[s] != SLOT_EMPTY)
				continue;
			if (refill(coll_q) == 0)
				break;                     /* nothing in hand and nothing queued */
			CollisionCandidate c = pending[first++];
			n_pending -= 1;
			assert(c.i == shared.header.i);
			seed0[s] = c.seed0;
			seed1[s] = c.seed1;
			len0[s] = c.len0_maybe;            /* 0 == it saturated: MEASURE counts it up from there */
			len1[s] = c.len1_maybe;
			measuring0[s] = (c.len0_maybe == 0);
			measuring1[s] = (c.len1_maybe == 0);
			end[s] = c.end;
			start0[s] = (shared.header.root_seed + params.multiplier * c.seed0) & wrapper.out_mask;
			start1[s] = (shared.header.root_seed + params.multiplier * c.seed1) & wrapper.out_mask;
			assert(not is_distinguished_point(start0[s], params.threshold));
			assert(not is_distinguished_point(start1[s], params.threshold));
			lane0[s] = free_lane[--n_free];
			lane1[s] = free_lane[--n_free];
			x[lane0[s]] = start0[s];
			x[lane1[s]] = start1[s];
			n_busy += 1;
			started += 1;
			if (measuring0[s])
				ctr[N_MEASURE] += 1;
			if (measuring1[s])
				ctr[N_MEASURE] += 1;
			if (measuring0[s] || measuring1[s])
				phase[s] = SLOT_MEASURE;
			else
				begin_align(s, ctr);       /* it picks the chain that has to move, or marches at once */
		}
		return started;
	}

	/*
	 * One vmixf over every lane, then one step of every slot.  N_EVAL counts the evaluations a scalar
	 * resolver would have made, not the lanes burnt, so that it stays comparable across the two paths.
	 */
	void step(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], const Params &params,
	          SharedContext<Scheme> &shared)
	{
		wrapper.vmixf(shared.header.i, x, y);

		for (int s = 0; s < nslots; s++) {
			if (phase[s] == SLOT_MEASURE) {
				if (measuring0[s] && not measure_chain(s, lane0[s], len0[s], measuring0[s], params, ctr))
					continue;          /* the slot is gone: its lanes are back in the pool */
				if (measuring1[s] && not measure_chain(s, lane1[s], len1[s], measuring1[s], params, ctr))
					continue;
				if (not measuring0[s] && not measuring1[s])
					begin_align(s, ctr);
			} else if (phase[s] == SLOT_ALIGN) {
				int l = (len0[s] > len1[s]) ? lane0[s] : lane1[s];   /* the longer chain is the mover */
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
					retire_collision(wrapper, ctr, hll, shared, seed0[s], seed1[s],
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
 * Retire up to `budget` candidates (0 == no limit) from `coll_q`, the queue of this walker's own
 * inserter.  A step costs a whole vmixf however few slots are in flight, so in the steady state we
 * stop as soon as nothing more is in hand and the batch is not full, and let the rest wait for the
 * next chunk: running the batch down to its last candidate would spend full-width steps on one or two
 * live lanes.  `drain` is the end of the round, where every candidate must be retired whatever it
 * costs.  `trail` is the scalar path's recording buffer, and is unused when vlen > 1.
 * PROTOCOL.md §3.2, §3.3.
 */
template <class ProblemWrapper>
void service_collisions(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], const Params &params,
                        SharedContext<Scheme> &shared, CollisionQueue &coll_q,
                        VecResolver<ProblemWrapper> &resolver, size_t budget, bool drain, u64 trail[])
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
 * PCS's producer, the walker: walks vlen trails in lockstep, ships every DP to the comm thread over its SPSC
 * queue, and between chunks retires the candidates queued by the one inserter of its own thread group
 * (ctx.group, PROTOCOL.md §1).  Chain indices are strided by n_producers from the global walker index.
 * Wind-down: ctx.state, PROTOCOL.md §3.3.
 */
template <class ProblemWrapper>
void Scheme::producer_thread(ThreadContext<Scheme> &ctx, const ProblemWrapper &wrapper, const Params &params,
                             SharedContext<Scheme> &shared, int walker_index)
{
	constexpr int vlen = ProblemWrapper::vlen;
	SPSCQueue &out = *ctx.q;
	CollisionQueue &coll_q = *shared.scheme.coll_q[ctx.group];
	u64 *ctr = ctx.ctr;
	u8 *hll = ctx.scheme.hll.data();

	const u64 i = shared.header.i;
	const u64 root_seed = shared.header.root_seed;

	int jbits = params.jbits;
	u64 jmask = make_mask(jbits);

	/* state of the vlen chains being walked */
	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 len[vlen];                  /* steps walked so far */
	u64 seed[vlen];                 /* the chain index j of each */

	/* the candidates being resolved, vlen/2 of them at once (PROTOCOL.md §3.2) */
	VecResolver<ProblemWrapper> resolver;

	/* the trail a scalar resolution records: the longest one a round can produce, plus its start */
	std::vector<u64> trail;
	if constexpr (vlen == 1)
		trail.resize(params.dp_max_it + 1);

	u64 j = (u64) params.rank * params.producers_per_node + walker_index;
	for (int k = 0; k < vlen; k++)
		start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, params.n_producers, k);
	assert((j & jmask) == j);

	for (;;) {
		int st = ctx.state.load(std::memory_order_acquire);

		if (st == HOLD) {
			ctx.state.store(HELD, std::memory_order_release);
			continue;
		}

		if (st == DRAIN) {
			service_collisions(wrapper, ctr, hll, params, shared, coll_q, resolver, 0, true,
			                   trail.data());
			assert(resolver.n_busy == 0);
			assert(resolver.n_pending == 0);
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}

		service_collisions(wrapper, ctr, hll, params, shared, coll_q, resolver, params.coll_per_chunk,
		                   false, trail.data());

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
					u64 l = std::min(len[k], params.len_sat);
					Point p = {x[k], (seed[k] & jmask) | (l << jbits)};
					if (not out.push(p))
						ctr[DROP_PRODUCERQ] += 1;
				}
				if (dp || failure) {
					if (failure && not dp)
						ctr[BAD_DP] += 1;
					start_chain(params, wrapper.out_mask, root_seed, j,
					            x, len, seed, params.n_producers, k);
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
