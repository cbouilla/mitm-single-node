#ifndef MITM_COMMON
#define MITM_COMMON

#include <string>
#include <cmath>
#include <climits>

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
    double beta = 10;             /* use function variant for beta*w distinguished points */
    double theta = -1;            /* proportion of distinguished points. -1 == auto-choose */
    u64 multiplier = 0x2545f4914f6cdd1dull;

    /* other relevant quantities deduced from the above */
    u64 threshold;                /* any integer less than this is a DP */
    u64 w;                        /* # slots in the dict */
    u64 dp_max_it;                /* how many iterations to find a DP. */
    u64 points_per_version;       /* #DP per version of the function */

	/* utilities */
    bool verbose = 1;             /* print progress information */

    static double optimal_theta(double w, int n)
    {
        return 2.25 * std::sqrt((double) w / (1ll << n));
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
	u64 colliding_len = 0;          // since beginning
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
	u64 colliding_len_i = 0;

	/* display management */
	bool display_active = 1;
 	u64 w;                          // #slots
 	int pb_n;
	u64 n_dp_prev = 0;              // #DP found since last display
	u64 n_points_prev = 0;          // #function eval since last display
	u64 bytes_sent_prev = 0;
	double last_display;
	
	Counters() {}
	Counters(bool display_active) : display_active(display_active) {}

	void ready(int n, u64 _w)
	{
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

	void round_display()
	{
		printf("\n%.2f avg trail length (%.2f collliding).  %.2f%% probe failure.  %.2f%% walk-robinhhod.  %.2f%% walk-noncolliding.  %.2f%% same-value\n",
                (double) n_points_trails / n_dp, 
                (double) colliding_len / 2.0 / n_collisions, 
                100. * bad_probe / n_dp_i, 
                100. * bad_walk_robinhood / n_dp_i, 
                100. * bad_walk_noncolliding / n_dp_i, 
                100. * bad_collision / n_dp_i);
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

	void found_collision(u64 len0, u64 len1) 
	{
		colliding_len += len0 + len1;
		colliding_len_i += len0 + len1;
		n_collisions += 1;
		n_collisions_i += 1;
	}

	// call this when the dictionnary is flushed / a new mixing function tried
	void flush_dict()
	{
		last_display = 0;
		display();
		round_display();
		n_flush += 1;
		n_dp_i = n_collisions_i = colliding_len_i = 0;
		last_update = wtime();
		bad_dp = bad_probe = bad_walk_robinhood = bad_walk_noncolliding = bad_collision = 0;
	}

	void done()
	{
		end_time = wtime();
		if (not display_active)
			return;
		double total_time = end_time - start_time;
		printf("\n----------------------------------------\n");
		printf("Took %0.2f sec to find the golden inputs\n", total_time);
		printf("Used %" PRId64 " ≈ 2^%0.2f mixing functions\n", 1+n_flush, std::log2(1+n_flush));
		// printf("Evaluated f() %" PRId64 " ≈ 2^%0.2f times\n", n_points, std::log2(n_points));
		// printf("  - %" PRId64 " ≈ 2^%0.2f times to find DPs (%.1f%%)\n", 
		// 	n_points_trails, std::log2(n_points_trails), 100.0 * n_points_trails / n_points);
		// u64 n_points_walk = n_points - n_points_trails;
		// printf("  - %" PRId64 " ≈ 2^%0.2f times in walks (%.1f%%)\n", 
		// 	n_points_walk, std::log2(n_points_walk), 100.0 * n_points_walk / n_points);
		printf("Found %" PRId64 " ≈ 2^%0.2f distinguished points\n", n_dp, std::log2(n_dp));
		printf("Found %" PRId64 " ≈ 2^%0.2f collisions\n", n_collisions, std::log2(n_collisions));
  }
};  
}
#endif