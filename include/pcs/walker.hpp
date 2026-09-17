#ifndef MITM_PCS_WALKER
#define MITM_PCS_WALKER

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <utility>

#include "tools.hpp"
#include "router/router.hpp"
#include "pcs/params.hpp"
#include "pcs/shared.hpp"

/*
 * PCS's producer: one machine of vlen lanes that walks trails and resolves collisions on the same lanes.
 * One vmixf per cycle advances every lane; a lane runs one instruction to its end -- a fresh trail to its
 * distinguished point, or the resolution of one candidate its dict thread queued -- and a lane that ends
 * one is redispatched at once, a pending collision before a fresh trail.  Nothing is ever preempted.
 *
 * Why one machine: a round spends a fair share of its evaluations locating collisions, and a resolver
 * that borrows the lanes in batches leaves a third of them idle (only the march moves both chains) and
 * steps its last candidates at full width.  Here a lane a collision does not need walks a trail instead.
 * A collision holds ONE lane for its whole life: its march alternates the two chains on it, which costs
 * exactly what two lanes would and leaves nothing to allocate or wait for.
 */

namespace mitm::pcs {

inline bool is_distinguished_point(u64 x, u64 threshold)
{
	return x <= threshold;
}

/*
 * A located collision (x0, x1), inputs of the same evaluation: tally it, fold it into this walker's
 * HyperLogLog, test the pair.  True == it was the golden pair, already published.
 */
template<class ProblemWrapper>
bool retire_collision(const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], Shared &shared,
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


/********************************* the lanes **********************************/

/*
 * A walker's lanes and what each is doing.  A trail walks from its start to a distinguished point.  A
 * collision walks its two chains through, on one lane:
 *
 *   MEASURE0/1  a length saturated on the wire or in the shard: re-walk that chain to its endpoint to learn it
 *   ALIGN       rewind both chains and step the longer one until both are the same distance from the endpoint
 *   MARCH       step both, chain 0 on the even cycles and chain 1 on the odd ones, and compare the images
 *               until they meet or the shorter trail runs out
 *
 * Every phase advances one chain by one evaluation per cycle, so a live lane never wastes an evaluation and
 * N_EVAL is the count of live lanes; a lane sits idle only in the drain, once nothing is left to start.
 * Every position is an output of the function or a start point, so a free lane's 0 keeps vmixf in range.
 */
template<class ProblemWrapper>
struct alignas(sizeof(u64) * ProblemWrapper::vlen) Lanes {
	static constexpr int vlen = ProblemWrapper::vlen;
	enum lane_phase {FREE, TRAIL, MEASURE0, MEASURE1, ALIGN, MARCH};

	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));   /* the point each lane evaluates this cycle */
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));   /* one vmixf of x */
	int phase[vlen];          /* what each lane is doing */
	u64 seed0[vlen];          /* chain 0's index: a trail's own, or the candidate's incoming point */
	u64 seed1[vlen];          /* chain 1's index: the point that was in the slot */
	u64 start0[vlen];         /* where chain 0 starts */
	u64 start1[vlen];         /* where chain 1 starts */
	u64 len0[vlen];           /* chain 0's trail length: a trail's steps so far; 0 == unknown, counted up in MEASURE0 */
	u64 len1[vlen];           /* chain 1's, likewise, counted up in MEASURE1 */
	u64 end[vlen];            /* the shared endpoint, in full: what a measured chain must reach */
	u64 remaining[vlen];      /* steps left in ALIGN, then in MARCH */
	u64 other[vlen];          /* MARCH: the position of the chain not on the lane this cycle */
	u64 y0[vlen];             /* MARCH: chain 0's image, saved from the even cycle */
	int turn[vlen];           /* MARCH parity: 0 == chain 0 is on the lane */
	int n_live;               /* lanes not FREE: what vmixf works for */

	/* the round's context, the walker's own */
	Router_thread &rt;              /* where the distinguished points go */
	const ProblemWrapper &wrapper;  /* the function the lanes walk */
	u64 *ctr;                       /* this walker's tallies */
	u8 *hll;                        /* this walker's distinct-collision registers */
	const Params &params;           /* the difficulty, the lengths' widths, the chain-index multiplier */
	Shared &shared;                 /* the round's header and the golden pair */
	const u64 jmask;                /* the low jbits: what a chain index ships with a point */
	const u64 n_recv;               /* the routing modulus: dict threads over all nodes */

	Lanes(Router_thread &rt, const ProblemWrapper &wrapper, u64 ctr[], u8 hll[], const Params &params, Shared &shared)
		: rt(rt), wrapper(wrapper), ctr(ctr), hll(hll), params(params), shared(shared),
		  jmask(make_mask(params.jbits)), n_recv((u64) Router_size(rt, ROUTER_RECEIVER, ROUTER_GLOBAL))
	{
		for (int l = 0; l < vlen; l++) {
			x[l] = 0;
			phase[l] = FREE;
			seed0[l] = 0;
			seed1[l] = 0;
			start0[l] = 0;
			start1[l] = 0;
			len0[l] = 0;
			len1[l] = 0;
			end[l] = 0;
			remaining[l] = 0;
			other[l] = 0;
			y0[l] = 0;
			turn[l] = 0;
		}
		n_live = 0;
	}

	/* the lane is done with its instruction, finished or abandoned.  Every such path ends here */
	void release(int l)
	{
		x[l] = 0;
		phase[l] = FREE;
		n_live -= 1;
	}

	/* a fresh trail on lane l: chain j, from `start`, which is not a distinguished point */
	void begin_trail(int l, u64 j, u64 start)
	{
		assert(not is_distinguished_point(start, params.threshold));
		x[l] = start;
		len0[l] = 0;
		seed0[l] = j;
		phase[l] = TRAIL;
		n_live += 1;
	}

	/*
	 * Both chains the same distance from their endpoint: step them together, chain 0 on the lane first.
	 * Equal positions here mean one trail is a suffix of the other, which is no collision.
	 */
	void begin_march(int l)
	{
		if (len0[l] < len1[l]) {          /* chain 1 was the mover: chain 0, still at its start, takes the lane */
			other[l] = x[l];
			x[l] = start0[l];
		} else {
			other[l] = start1[l];
		}
		if (x[l] == other[l]) {
			ctr[BAD_WALK_ROBINHOOD] += 1;
			release(l);
			return;
		}
		remaining[l] = std::min(len0[l], len1[l]);
		turn[l] = 0;
		phase[l] = MARCH;
	}

	/* both lengths exact: rewind, and step the longer chain down to the shorter one's length */
	void begin_align(int l)
	{
		assert(len0[l] > 0 && len1[l] > 0);
		if (len0[l] >= len1[l]) {
			x[l] = start0[l];
			remaining[l] = len0[l] - len1[l];
		} else {
			x[l] = start1[l];
			remaining[l] = len1[l] - len0[l];
		}
		phase[l] = ALIGN;
		if (remaining[l] == 0)
			begin_march(l);
	}

	/* chain 0's length is exact: measure chain 1 if its own is not, else align */
	void after_chain0(int l)
	{
		if (len1[l] == 0) {
			ctr[N_MEASURE] += 1;
			x[l] = start1[l];
			phase[l] = MEASURE1;
		} else {
			begin_align(l);
		}
	}

	/* a candidate on lane l, at the first phase it needs */
	void begin_collision(int l, const CollisionCandidate &c)
	{
		const u64 root_seed = shared.header.root_seed;
		assert(c.i == shared.header.i);
		seed0[l] = c.seed0;
		seed1[l] = c.seed1;
		len0[l] = c.len0_maybe;            /* 0 == it saturated: MEASURE counts it up from there */
		len1[l] = c.len1_maybe;
		end[l] = c.end;
		start0[l] = (root_seed + params.multiplier * c.seed0) & wrapper.out_mask;
		start1[l] = (root_seed + params.multiplier * c.seed1) & wrapper.out_mask;
		assert(not is_distinguished_point(start0[l], params.threshold));
		assert(not is_distinguished_point(start1[l], params.threshold));
		n_live += 1;
		if (len0[l] == 0) {
			ctr[N_MEASURE] += 1;
			x[l] = start0[l];
			phase[l] = MEASURE0;
		} else {
			after_chain0(l);
		}
	}

	/*
	 * One step of a chain re-walked to its distinguished point to recover its length.  True once it is
	 * there and the length is exact; false while it walks, or when the candidate is gone and the lane free.
	 */
	bool measure_step(int l, u64 &len)
	{
		len += 1;
		x[l] = y[l];
		if (is_distinguished_point(x[l], params.threshold)) {
			if (x[l] == end[l])
				return true;
			ctr[BAD_WALK_NONCOLLIDING] += 1;      /* not the trail the dictionary meant */
			release(l);
			return false;
		}
		if (len == params.dp_max_it) {
			ctr[BAD_DP] += 1;                     /* the re-walk never reached a distinguished point */
			release(l);
		}
		return false;
	}

	/* one vmixf over every lane, then what each live lane does with its image */
	void step()
	{
		wrapper.vmixf(shared.header.i, x, y);
		ctr[N_EVAL] += n_live;

		for (int l = 0; l < vlen; l++) {
			switch (phase[l]) {
			case FREE:
				break;
			case TRAIL:
				len0[l] += 1;
				x[l] = y[l];
				if (is_distinguished_point(x[l], params.threshold)) {
					ctr[N_DP] += 1;
					ctr[N_POINTS_TRAILS] += len0[l];
					Router_Push(x[l], (seed0[l] & jmask) | (std::min(len0[l], params.len_sat) << params.jbits),
					            (int) (x[l] % n_recv), rt);
					release(l);
				} else if (len0[l] == params.dp_max_it) {
					ctr[BAD_DP] += 1;
					release(l);
				}
				break;
			case MEASURE0:
				if (measure_step(l, len0[l]))
					after_chain0(l);
				break;
			case MEASURE1:
				if (measure_step(l, len1[l]))
					begin_align(l);
				break;
			case ALIGN:
				x[l] = y[l];
				remaining[l] -= 1;
				if (remaining[l] == 0)
					begin_march(l);
				break;
			case MARCH:
				if (turn[l] == 0) {                    /* even cycle: chain 0's image; chain 1 takes the lane */
					y0[l] = y[l];
					std::swap(x[l], other[l]);
					turn[l] = 1;
				} else if (y[l] == y0[l]) {            /* odd cycle: the two images meet */
					retire_collision(wrapper, ctr, hll, shared, seed0[l], seed1[l], other[l], x[l], len0[l], len1[l]);
					release(l);
				} else {                               /* odd cycle: both advance, chain 0 back on the lane */
					other[l] = y[l];
					x[l] = y0[l];
					turn[l] = 0;
					remaining[l] -= 1;
					if (remaining[l] == 0) {
						ctr[BAD_WALK_NONCOLLIDING] += 1;   /* the trails never met: a false positive */
						release(l);
					}
				}
				break;
			}
		}
	}
};


/********************************* the round **********************************/

/*
 * A walker's round.  Every cycle: a free lane takes a pending collision if there is one, a fresh trail
 * otherwise; then one step of the machine.  Candidates come off the queue in runs of up to 64, so the lock
 * is taken once per run; chain indices are strided by the number of walkers from this one's global
 * index, so no two walkers ever share a chain and the partition needs no message.
 *
 * The controller ends the round: `round_over` closes this walker and abandons its trails, and what
 * follows is the drain, where every candidate left must be resolved before the next round mixes
 * differently.  A walker may not leave before its dict thread has stopped pushing AND the queue it
 * shares with its fellows is empty -- `done` first, then the queue, or a run pushed between the two
 * would be left behind.
 */
template <class ProblemWrapper>
void walker_round(Router_thread &rt, const ProblemWrapper &wrapper, const Params &params, Shared &shared,
                  u64 *ctr, u8 *hll, int r)
{
	constexpr int vlen = ProblemWrapper::vlen;
	constexpr int PENDING = 64;                             /* candidates taken off the queue at once */
	CollisionQueue &coll_q = *shared.chan[r].q;
	const u64 root_seed = shared.header.root_seed;
	Lanes<ProblemWrapper> lanes(rt, wrapper, ctr, hll, params, shared);   /* every lane free; the round's own */
	CollisionCandidate pending[PENDING];                    /* the run this walker took from the queue */
	int n_pending = 0;                                      /* how many are left in the run */
	u64 j = (u64) Router_rank(rt, ROUTER_GLOBAL);           /* the chain-index cursor: this walker's rank among all walkers */
	const u64 jinc = (u64) Router_size(rt, ROUTER_SENDER, ROUTER_GLOBAL);   /* its stride */
	bool closed = false;

	for (;;) {
		if (not closed && shared.round_over) {
			Router_Close(rt);
			for (int l = 0; l < vlen; l++)             /* what a trail in flight has walked is lost: its lane goes to the collisions */
				if (lanes.phase[l] == Lanes<ProblemWrapper>::TRAIL)
					lanes.release(l);
			closed = true;
		}
		if (lanes.n_live < vlen) {
			for (int l = 0; l < vlen; l++) {
				if (lanes.phase[l] != Lanes<ProblemWrapper>::FREE)
					continue;
				if (n_pending == 0 && not coll_q.is_empty())     /* relaxed read: an idle queue costs no lock traffic */
					n_pending = (int) coll_q.pop_bulk(pending, PENDING);
				if (n_pending > 0) {
					n_pending -= 1;
					lanes.begin_collision(l, pending[n_pending]);
				} else if (not closed) {
					u64 start;                                   /* the next chain of this walker's that is not itself a distinguished point */
					for (;;) {
						j += jinc;
						start = (root_seed + j * params.multiplier) & wrapper.out_mask;
						if (not is_distinguished_point(start, params.threshold))
							break;
					}
					assert((j & lanes.jmask) == j);              /* the chain index must fit beside the length */
					lanes.begin_trail(l, j, start);
				}
			}
		}
		if (lanes.n_live == 0) {
			if (closed && shared.chan[r].done.load_acquire() && coll_q.is_empty())
				break;
			cpu_relax();
			continue;
		}
		lanes.step();
	}
	assert(lanes.n_live == 0);
	assert(n_pending == 0);
}

}
#endif
