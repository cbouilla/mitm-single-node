#ifndef MITM_COUNTERS
#define MITM_COUNTERS

#include <cmath>
#include <cassert>
#include <strings.h>          // ffsll, for the HyperLogLog

#include "tools.hpp"

namespace mitm {

/*
 * Every diagnostic tally of the engine, by index.  One array of u64 holds them all,
 * whichever thread produces them, so that the same layout serves for the per-thread
 * tallies, the progress reports (deltas, TAG_REPORT) and the end-of-round MPI_Reduce.
 *
 * None of this is part of the attack: it is what tells us whether the parameters are
 * any good.
 */
enum counter {
	/* walkers */
	N_EVAL = 0,             /* evaluations of the mixing function, by walking or by resolving */
	N_DP,                   /* distinguished points found */
	N_POINTS_TRAILS,        /* sum of the lengths of the trails that reached a DP */
	N_COLLISIONS,           /* collisions located */
	COLLIDING_LEN_MIN,      /* sum of the shorter length of each colliding pair */
	COLLIDING_LEN_MAX,      /* ... and of the longer one */
	BAD_DP,                 /* trail gave up before reaching a distinguished point */
	BAD_COLLISION,          /* the two trails "collide" on the same value */
	BAD_WALK_ROBINHOOD,     /* one trail is a suffix of the other */
	BAD_WALK_NONCOLLIDING,  /* dictionary false positive: the trails never meet */
	DROP_WALKERQ,           /* DP dropped: the queue to the comm thread was full */
	/* inserters */
	N_PROBE,                /* dictionary probes retired */
	BAD_PROBE,              /* dictionary slot was empty, or held a different key */
	DROP_COLL,              /* candidate dropped: the collision queue was full */
	/* comm thread */
	DROP_OUT,               /* DP dropped: the outgoing MPI buffer was still in flight */
	DROP_INSERTERQ,         /* DP dropped: a local inserter's queue was full */
	N_COUNTERS
};

/*
 * The tallies of ONE round (= one version of the mixing function).  Every thread owns
 * one and adds to it with no synchronisation at all: plain u64, not atomics.  The
 * comm thread reads them twice: during the round, after an `omp flush`, to build its
 * progress reports -- a count that lags by a unit or two is fine there -- and after
 * the end-of-round barrier, exactly, to merge them, reduce the result across the
 * nodes and hand it to the controller.  All-time totals and every printf live there.
 */
class Counters {
public:
	u64 c[N_COUNTERS] = {};

	/* HyperLogLog over the collisions of this round, used to estimate how many of
	   them are DISTINCT -- the quantity that actually drives the attack. */
	vector<u8> hll;

	Counters() : hll(0x10000) {}

	// collect statistics
	void found_collision(u64 x0, u64 len0, u64 x1, u64 len1)
	{
		if (len0 < len1) {
			c[COLLIDING_LEN_MIN] += len0;
			c[COLLIDING_LEN_MAX] += len1;
		} else {
			c[COLLIDING_LEN_MIN] += len1;
			c[COLLIDING_LEN_MAX] += len0;
		}
		c[N_COLLISIONS] += 1;

		u64 h = murmur128(x0, x1);
		u64 idx = h >> 48;
		assert(idx < 0x10000);
		int rho = ffsll(h);
		if (hll[idx] < rho)
			hll[idx] = rho;
	}

	/*
	 * Fold another set of tallies into this one.  Sums the counts and takes the
	 * element-wise max of the HyperLogLog registers, which is exactly the HLL merge
	 * rule -- the same operations as the MPI_SUM / MPI_MAX reductions across nodes,
	 * and the controller folds each round into its all-time totals the same way.
	 */
	void merge(const Counters &o)
	{
		for (int k = 0; k < N_COUNTERS; k++)
			c[k] += o.c[k];
		for (int k = 0; k < 0x10000; k++)
			if (hll[k] < o.hll[k])
				hll[k] = o.hll[k];
	}

	void reset()
	{
		for (int k = 0; k < N_COUNTERS; k++)
			c[k] = 0;
		hll.assign(0x10000, 0);
	}

	/* uses the HyperLogLog algorithm */
	static u64 distinct_collisions_estimation(const vector<u8> &h)
	{
		double acc = 0;
		double alpha = 0.7213 / (1 + 1.079 / 0x10000);
		for (int i = 0; i < 0x10000; i++)
			acc += 1.0 / (1 << h[i]);
		double E = alpha * 0x100000000 / acc;
		if (E >= 2.5 * 0x10000)
			return E;
		// low cardinality, potential correction
		int V = 0;
		for (int i = 0; i < 0x10000; i++)
			if (h[i] == 0)
				V += 1;
		if (V == 0)
			return E;
		else
			return 0x10000 * log(65536.0 / V);
	}
};

}
#endif
