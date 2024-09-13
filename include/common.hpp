#ifndef MITM_COMMON
#define MITM_COMMON

#include <string>
#include <cmath>
#include <climits>
#include <cstring>

// base classes for PCS and naive algorithm

#include "tools.hpp"
#include "dict.hpp"

namespace mitm {

class Parameters {
public:
    /* hardware-dependent */
    u64 nbytes_memory = 0;        /* how much RAM to use on each machine */
    int n_nodes = 1;              /* #hosts (with shared RAM) */
    int n_recv = 1;               /* #instances of the dictionary */

    /* algorithm parameters */
    double alpha = 2.5;           /* auto-chosen theta == alpha * sqrt(w/n) */
    double beta = 8;              /* use function variant for beta*w distinguished points */
    double theta = -1;            /* proportion of distinguished points. -1 == auto-choose */

    u64 multiplier = 0x2545f4914f6cdd1dull;       /* to generate starting points */

    /* other relevant quantities deduced from the above */
    u64 threshold;                /* any integer less than this is a DP */
    u64 w;                        /* # slots in the dict */
    u64 dp_max_it;                /* how many iterations to find a DP. */
    u64 points_per_version;       /* #DP per version of the function */

	/* utilities */
    bool verbose = 1;             /* print progress information */
    u64 max_versions = 0xffffffffffffffffull;       /* how many functions to try before giving up */


    double optimal_theta(double w, int n)
    {
        return alpha * std::sqrt((double) w / (1ll << n));
    }

    void finalize(int n, int m)
    {
        if (nbytes_memory == 0)
            errx(1, "the amount of RAM to use (per node) must be specified");

        w = PcsDict::get_nslots(nbytes_memory * n_nodes, n_recv);
        /* auto-choose the difficulty if not set */
        double auto_theta = optimal_theta(w, n);
        if (theta < 0) {
        	theta = auto_theta;
            if (theta > 1)
                theta = 1;
            if (verbose)
                printf("AUTO-TUNING: setting 1/theta == %.2f\n", 1 / theta);
        } else {
            if (verbose)
                printf("NOTICE: using 1/theta == %.2f vs ``optimal'' 1/theta == %.2f\n", 1/theta, auto_theta);
        }
        threshold = pow(2, m) * theta;
        dp_max_it = 20 / theta;
        points_per_version = beta * w;

        /* display warnings if problematic choices were made */
        if (verbose && theta == 1) {
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
            printf("---> zero difficulty (use the naive technique!)\n");
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");            
        }
    }
};


/* Non-essential counters but helpful to have, e.g. n_collisions/sec */
class Counters {
public:
 	const u64 interval = 1LL << 16; // for printing only
 	const double min_delay = 1;   // for printing only
 	
 	/* stats for the full computation */
	u64 n_flush = 0;                // #times dict were flushed
	u64 n_dp = 0;                   // since beginning
	u64 n_points_trails = 0;
	u64 n_collisions = 0;           // since beginning
	u64 colliding_len_min = 0;          // since beginning
	u64 colliding_len_max = 0;          // since beginning
	u64 bad_dp = 0;
	u64 bad_probe = 0;
	u64 bad_collision = 0;
	u64 bad_walk_robinhood = 0;
	u64 bad_walk_noncolliding = 0;
	double start_time;
	double end_time;

	/* stats for this i */
	double last_update;             // timestamp of the last flush
	u64 n_dp_i = 0;                     // since last flush
	u64 n_collisions_i = 0;             // since last flush
	u64 n_coll_unique = 0;             // since last flush
	u64 colliding_len_min_i = 0;
	u64 colliding_len_max_i = 0;

	/* display management */
	bool display_active = 1;
 	u64 w;                          // #slots
 	int pb_n;
	u64 n_dp_prev = 0;              // #DP found since last display
	u64 n_points_prev = 0;          // #function eval since last display
	u64 bytes_sent_prev = 0;
	double last_display;
	
	vector<u8> hll, hll_i;

	Counters() {}
	Counters(bool display_active) : display_active(display_active) {}

	void ready(int n, u64 _w)
	{
		hll.resize(0x10000);
		hll_i.resize(0x10000);
		pb_n = n;
		w = _w;
		double now = wtime();
		start_time = now;
		start_time = now;
		last_display = now;
		last_update = now;
	}

	void dp_failure()
	{ 
		bad_dp += 1; 
	}
	
	void probe_failure()
	{
		bad_probe += 1;
	}
	
	void walk_robinhood() {
		bad_walk_robinhood += 1;
	}

	void walk_noncolliding() {
		bad_walk_noncolliding += 1;
	}
	
	void collision_failure() {
		bad_collision += 1;
	}

	// call this when a new DP is found
	void found_distinguished_point(u64 chain_len)
	{
		n_dp += 1;
		n_dp_i += 1;
		n_points_trails += chain_len;
		if ((n_dp % interval == interval - 1))
			 display();
	}

	void found_collision(u64 x0, u64 len0, u64 x1, u64 len1) 
	{
		if (len0 < len1) {
			colliding_len_min += len0;
			colliding_len_max += len1;
			colliding_len_min_i += len0;
			colliding_len_max_i += len1;
		} else {
			colliding_len_min += len1;
			colliding_len_max += len0;
			colliding_len_min_i += len1;
			colliding_len_max_i += len0;
		}
		n_collisions += 1;
		n_collisions_i += 1;

		/* update HyperLogLog */
    	u64 h = murmur128(x0, x1);
    	u64 idx = h >> 48;
    	assert(idx < 0x10000);
    	int rho = ffsll(h);
    	if (hll[idx] < rho)
    	    hll[idx] = rho;
    	if (hll_i[idx] < rho)
    	    hll_i[idx] = rho;
	}

	/* uses the HyperLogLog algorithm */
	static u64 distinct_collisions_estimation(const vector<u8> h)
	{
		double acc;
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

	// call this when the dictionnary is flushed / a new mixing function tried
	void flush_dict()
	{
		last_display = 0;
		display();
		round_display();
		n_flush += 1;
		n_dp_i = n_collisions_i = colliding_len_min_i = colliding_len_max_i = 0;
		last_update = wtime();
		bad_dp = bad_probe = bad_walk_robinhood = bad_walk_noncolliding = bad_collision = 0;
		n_coll_unique += distinct_collisions_estimation(hll_i);
		hll_i.clear();
		hll_i.resize(0x10000);
	}

	/************************** verbosity ************************/

	// invoked regularly
	void display()
	{
		if (not display_active)
			return;
		double now = wtime();
		double delta = now - last_display;
		if (delta >= min_delay) {
			u64 N = 1ull << pb_n;
			double dprate = (n_dp - n_dp_prev) / delta;
			char hdprate[8];
			human_format(dprate, hdprate);
			printf("\r#i = %" PRId64 " (%.02f*n/w). %s DP/sec.  #DP (this i / total) %.02f*w / %.02f*n.  #coll (this i / total) %.02f*w / %0.2f*n",
			 		n_flush, (double) n_flush * w / N, 
			 		hdprate, 
			 		(double) n_dp_i / w, (double) n_dp / N,
			 		(double) n_collisions_i / w, (double) n_collisions / N);
			fflush(stdout);
			n_dp_prev = n_dp;
			last_display = wtime();
		}
	}

/* 
 X_i == #red balls in bin i
 Y_i == #blue balls in bin i

 #collisions == sum_i X_i Y_i

 E[#collisions] == 2^m E[X_i * Y_i] == 2^m E[X_i] * E[Y_i] == 2^m (2^n / 2^m) * (2^n / 2^m) == 2^(2n - m)

*/
	// each new version of the function
	void round_display()
	{
		u64 N = 1ull << pb_n;
		u64 E_i = distinct_collisions_estimation(hll_i);
		u64 E = distinct_collisions_estimation(hll);
		double E_exp = N * (1 - std::pow(1 - 1. / N, n_collisions));
		double avglen = (double) n_points_trails / n_dp;
		printf("\n%.2f avg trail length (x%.2f & x%.2f collliding).  %.2f%% probe failure.  %.2f%% walk-robinhhod.  %.2f%% walk-noncolliding.  %.2f%% same-value\n",
                avglen, 
                (double) colliding_len_min / n_collisions / avglen, 
                (double) colliding_len_max / n_collisions / avglen, 
                100. * bad_probe / n_dp_i, 
                100. * bad_walk_robinhood / n_dp_i, 
                100. * bad_walk_noncolliding / n_dp_i, 
                100. * bad_collision / n_dp_i);
		printf("#coll (this i / distinct / total / distinct / expected) %.02f*w / %.02f*w / %.02f*n / %.02f*n / %.02f*n\n", 
				(double) n_collisions_i / w,
                (double) E_i / w,
                (double) n_collisions / N,
                (double) E / N,
                (double) E_exp / N);
		printf("\n");
		fflush(stdout);
	}

	// last call
	void done()
	{
		if (not display_active)
			return;
		double total_time = wtime() - start_time;
		u64 N = 1ull << pb_n;
		u64 E = distinct_collisions_estimation(hll);
		u64 v = n_flush;
		printf("\n----------------------------------------\n");
		printf("Total running time %0.2fs\n", total_time);
		printf("Used %" PRId64 " = %.2f*n/w mixing functions\n", v, (double) v / N * w);
		// printf("Evaluated f() %" PRId64 " ≈ 2^%0.2f times\n", n_points, std::log2(n_points));
		// printf("  - %" PRId64 " ≈ 2^%0.2f times to find DPs (%.1f%%)\n", 
		// 	n_points_trails, std::log2(n_points_trails), 100.0 * n_points_trails / n_points);
		// u64 n_points_walk = n_points - n_points_trails;
		// printf("  - %" PRId64 " ≈ 2^%0.2f times in walks (%.1f%%)\n", 
		// 	n_points_walk, std::log2(n_points_walk), 100.0 * n_points_walk / n_points);
		// printf("Found %" PRId64 " ≈ 2^%0.2f distinguished points\n", n_dp, std::log2(n_dp));
		printf("Found %.02f*n collisions (%.02f*n distinct)\n", (double) n_collisions / N, (double) E / N);
		printf("Found %.02f*w collisions / version (%.02f*w distinct)\n", (double) n_collisions / v / w, (double) n_coll_unique / v / w);
	}
};  
}
#endif