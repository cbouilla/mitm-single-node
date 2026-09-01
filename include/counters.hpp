#ifndef MITM_COUNTERS
#define MITM_COUNTERS

#include <cmath>
#include <cassert>
#include <strings.h>          // ffsll, for the HyperLogLog

#include "tools.hpp"

namespace mitm {

/*
 * Diagnostic tallies for ONE round (= one version of the mixing function).  Every
 * thread owns one and adds to it without any synchronisation; at the end of the round
 * the comm thread merges them, reduces the result across the nodes and hands it to the
 * controller.  All-time totals and every printf live there, so nothing here is ever
 * read while the round is running.
 *
 * None of this is part of the attack: it is what tells us whether the parameters are
 * any good.
 */
class Counters {
public:
	u64 n_eval = 0;                 // evaluations of the mixing function, by walking or by resolving
	u64 n_points_trails = 0;        // sum of the lengths of the trails that reached a DP
	u64 n_collisions = 0;
	u64 colliding_len_min = 0;      // sum of the shorter length of each colliding pair
	u64 colliding_len_max = 0;      // ... and of the longer one
	u64 bad_dp = 0;                 // trail gave up before reaching a distinguished point
	u64 bad_probe = 0;              // dictionary slot was empty, or held a different key
	u64 bad_collision = 0;          // the two trails "collide" on the same value
	u64 bad_walk_robinhood = 0;     // one trail is a suffix of the other
	u64 bad_walk_noncolliding = 0;  // dictionary false positive: the trails never meet

	/* HyperLogLog over the collisions of this round, used to estimate how many of
	   them are DISTINCT -- the quantity that actually drives the attack. */
	vector<u8> hll;

	Counters() : hll(0x10000) {}

	void found_collision(u64 x0, u64 len0, u64 x1, u64 len1)
	{
		if (len0 < len1) {
			colliding_len_min += len0;
			colliding_len_max += len1;
		} else {
			colliding_len_min += len1;
			colliding_len_max += len0;
		}
		n_collisions += 1;

		u64 h = murmur128(x0, x1);
		u64 idx = h >> 48;
		assert(idx < 0x10000);
		int rho = ffsll(h);
		if (hll[idx] < rho)
			hll[idx] = rho;
	}

	/*
	 * Fold another thread's tallies into this one.  Sums the counts and takes the
	 * element-wise max of the HyperLogLog registers, which is exactly the HLL merge
	 * rule -- and the same operation as the MPI_MAX reduction used across nodes.
	 */
	void merge(const Counters &o)
	{
		n_eval += o.n_eval;
		n_points_trails += o.n_points_trails;
		n_collisions += o.n_collisions;
		colliding_len_min += o.colliding_len_min;
		colliding_len_max += o.colliding_len_max;
		bad_dp += o.bad_dp;
		bad_probe += o.bad_probe;
		bad_collision += o.bad_collision;
		bad_walk_robinhood += o.bad_walk_robinhood;
		bad_walk_noncolliding += o.bad_walk_noncolliding;

		for (int k = 0; k < 0x10000; k++)
			if (hll[k] < o.hll[k])
				hll[k] = o.hll[k];
	}

	void reset()
	{
		n_eval = n_points_trails = n_collisions = colliding_len_min = colliding_len_max = 0;
		bad_dp = bad_probe = bad_collision = bad_walk_robinhood = bad_walk_noncolliding = 0;
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
