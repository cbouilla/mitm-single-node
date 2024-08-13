#ifndef MITM_COUNTERS
#define MITM_COUNTERS

#include "common.hpp"

#include <string>
#include <vector>
#include <cmath>

namespace mitm {

/* Non-essential counters but helpful to have, e.g. n_collisions/sec */
struct Counters {
 	const u64 interval = 1LL << 16; // for printing only
 	const double min_delay = 1;   // for printing only
 	
 	/* stats for the full computation */
	u64 n_flush = 0;                // #times dict were flushed
	u64 n_points = 0;               // #function evaluation in total
	u64 n_points_trails = 0;        // #function evaluation to find DP
	u64 n_dp = 0;                   // since beginning
	u64 n_collisions = 0;           // since beginning
	u64 bad_dp = 0;
	u64 bad_probe = 0;
	u64 bad_walk = 0;
	u64 bad_collision = 0;
	double start_time;
	double end_time;

	/* waiting times for MPI implementation */
	double send_wait = 0.;          // time sender processes wait for the receiver to be ready
	double recv_wait = 0.;          // time the receivers wait for data from the senders
	u64 bytes_sent = 0;

	/* stats for this i */
	double last_update;             // timestamp of the last flush
	u64 n_dp_i = 0;                 // since last flush
	u64 n_collisions_i = 0;         // since last flush

	/* display management */
 	u64 w;                          // #slots
 	int pb_n;
	bool display_active = 1;
	u64 n_dp_prev = 0;              // #DP found since last display
	u64 n_points_prev = 0;          // #function eval since last display
	u64 bytes_sent_prev = 0;
	double last_display;
	
	Counters() {}
	Counters(bool display_active) : display_active(display_active) {}

	// only useful for MPI engine
	void reset()
	{
		n_points = 0;
		n_dp = 0;
		n_collisions = 0;
		send_wait = 0;
		recv_wait = 0;
		bytes_sent = 0;
	}

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
		#pragma omp atomic
		bad_dp += 1; 
	}
	
	void probe_failure()
	{
		#pragma omp atomic		
		bad_probe += 1;
	}
	
	void walk_failure() {
		#pragma omp atomic		
		bad_walk += 1;
	}
	
	void collision_failure() {
		#pragma omp atomic		
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
			double frate = (n_points - n_points_prev) / delta;
			char hfrate[8];
			human_format(frate, hfrate);
			double nrate = (bytes_sent - bytes_sent_prev) / delta;
			char hnrate[8];
			human_format(nrate, hnrate);
			printf("\r#i = %" PRId64 " (%.02f*n/w). %s f/sec.  %s DP/sec.  #DP (this i / total) %.02f*w / %.02f*n.  #coll (this i / total) %.02f*w / %0.2f*n.  Total #f=2^%.02f.  S/R wait %.02fs / %.02fs.  S->R %sB/s",
			 		n_flush, (double) n_flush * w / N, 
			 		hfrate, hdprate, 
			 		(double) n_dp_i / w, (double) n_dp / N,
			 		(double) n_collisions_i / w, (double) n_collisions / N,
			 		std::log2(n_points),
			 		send_wait, recv_wait,
			 		hnrate);
			fflush(stdout);
			n_dp_prev = n_dp;
			n_points_prev = n_points;
			bytes_sent_prev = bytes_sent;
			last_display = wtime();
		}
	}

	void function_evaluation()
	{
		#pragma omp atomic
		n_points += 1;
	}

	// call this when a new DP is found
	void found_distinguished_point(u64 chain_len)
	{
		#pragma omp atomic
		n_dp += 1;
		#pragma omp atomic
		n_dp_i += 1;
		#pragma omp atomic
		n_points_trails += chain_len;
		if ((n_dp % interval == interval - 1))
			 #pragma omp critical(counters)
			 display();
	}

	void found_collision() {
		#pragma omp atomic
		n_collisions += 1;
		#pragma omp atomic
		n_collisions_i += 1;
	}

	// call this when the dictionnary is flushed / a new mixing function tried
	void flush_dict()
	{
		display();
		/*double elapsed = wtime() - last_update;
		printf("\nUpdating the mixing function. The previous version:\n");
		printf(" - lasted %0.2f sec\n", elapsed);
		printf(" - 2^%0.2f function evaluations\n", std::log2(n_points));
		printf(" - generated %" PRId64 " ≈ 2^%0.2f distinguished points\n", n_dp[n_flush], std::log2(n_dp[n_flush]));
		printf(" - found %" PRId64 " ≈ 2^%0.2f collisions\n", n_collisions[n_flush], std::log2(n_collisions[n_flush]));
		printf(" - %" PRId64 " DP failure, %" PRId64 " probe failure, %" PRId64 " walk failure, %" PRId64 " collision failure\n", 
			bad_dp, bad_probe, bad_walk, bad_collision);
		printf("Currently totaling 2^%.02f iterations\n", std::log2(n_points));
		*/
		n_flush += 1;
		n_dp_i = 0;
		n_collisions_i = 0;
		last_update = wtime();
		bad_dp = bad_probe = bad_walk = bad_collision = 0;
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
		printf("Evaluated f() %" PRId64 " ≈ 2^%0.2f times\n", n_points, std::log2(n_points));
		printf("  - %" PRId64 " ≈ 2^%0.2f times to find DPs (%.1f%%)\n", 
			n_points_trails, std::log2(n_points_trails), 100.0 * n_points_trails / n_points);
		u64 n_points_walk = n_points - n_points_trails;
		printf("  - %" PRId64 " ≈ 2^%0.2f times in walks (%.1f%%)\n", 
			n_points_walk, std::log2(n_points_walk), 100.0 * n_points_walk / n_points);
		printf("Found %" PRId64 " ≈ 2^%0.2f distinguished points\n", n_dp, std::log2(n_dp));
		printf("Found %" PRId64 " ≈ 2^%0.2f collisions\n", n_collisions, std::log2(n_collisions));
  }
};  
}
#endif