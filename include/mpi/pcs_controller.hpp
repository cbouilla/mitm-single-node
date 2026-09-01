#ifndef MITM_MPI_CONTROLLER
#define MITM_MPI_CONTROLLER

#include <cmath>
#include <mpi.h>

#include "../common.hpp"
#include "../engine_common.hpp"
#include "common.hpp"
#include "pcs_comm.hpp"

namespace mitm {

/*
 * The controller lives on the comm thread of rank 0, which also has a full node to
 * run.  It is therefore a state object driven by that thread's poll loop rather than
 * a loop of its own: nothing here blocks, and nothing probes.
 */
class Controller {
public:
	const MpiParameters &params;

	u64 nround = 0;
	u64 stop = 0;
	optional<tuple<u64,u64,u64>> solution;         /* (i, x0, x1) */

	/* this round */
	u64 ndp = 0;                                   /* DPs reported by all nodes */
	int n_active_nodes = 0;
	double round_start = 0, last_display = 0;
	u64 drop_walkerq = 0, drop_out = 0, drop_inserterq = 0, drop_coll = 0;
	u64 nprobe = 0;                                /* dictionary probes retired */

	/* totals */
	u64 ndp_total = 0, ncoll_total = 0, nf_total = 0;
	double start_time;
	vector<u8> hll;                                /* all collisions since the start */

	Controller(const MpiParameters &p) : params(p), n_active_nodes(p.n_nodes), hll(0x10000)
	{
		start_time = wtime();
	}

	void banner(const PRNG &prng, u64 w_slots)
	{
		char hbuf[8], hdict[8];
		u64 bufbytes = (u64) 2 * params.n_nodes * DP_WORDS * sizeof(u64) * params.buffer_capacity;
		human_format(bufbytes, hbuf);
		human_format(params.n_nodes * w_slots * sizeof(u64), hdict);
		printf("Starting MPI+OpenMP collision search with seed=%016" PRIx64 "\n", prng.seed);
		printf("RAM per node == %sB buffers + dict.  Total dict == %sB (2^%.2f slots)\n",
			hbuf, hdict, std::log2((double) params.w));
		printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished points / version\n",
			params.beta, params.points_per_version, std::log2(params.points_per_version));
		fflush(stdout);
	}

	void begin_round()
	{
		ndp = 0;
		n_active_nodes = params.n_nodes;
		drop_walkerq = drop_out = drop_inserterq = drop_coll = nprobe = 0;
		round_start = wtime();
		last_display = round_start;
	}


	/*
	 * Digest one node report.  Returns the assignment to send back, or -1 when no
	 * reply is due (solution reports are one-way).
	 */
	int handle_report(const u64 r[REP_NWORDS])
	{
		if (r[REP_GOLDEN]) {
			if (not solution)
				solution = optional(tuple(r[REP_I], r[REP_X0], r[REP_X1]));
			stop = 1;
			return -1;
		}

		ndp += r[REP_NDP];
		drop_walkerq += r[REP_DROP_WALKERQ];
		drop_out   += r[REP_DROP_OUT];
		drop_inserterq += r[REP_DROP_INSERTERQ];
		drop_coll  += r[REP_DROP_COLL];
		nprobe     += r[REP_NPROBE];

		int assignment = KEEP_GOING;
		if (stop || ndp >= params.points_per_version) {
			assignment = NEW_VERSION;
			n_active_nodes -= 1;
		}
		display();
		return assignment;
	}

	void display()
	{
		double now = wtime();
		if (now - last_display <= 0.5)
			return;
		last_display = now;
		double delta = now - round_start;
		double dp_rate = ndp / delta;
		double completion = (double) ndp / params.points_per_version;
		char hrate[8], hnrate[8], hprobe[8];
		human_format(dp_rate / params.theta / params.n_walkers, hrate);
		human_format(ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);
		human_format((double) nprobe / params.n_inserters / delta, hprobe);
		printf("\rRound %" PRId64 ":  %.1fs (%.1f%%, ETA %.1fs).  %.2f*w #DP.  %s #f/s per walker.  %s probe/s per inserter.  node-->%sB/s   ",
			nround, delta, 100. * completion,
			(completion > 0) ? delta * (1 - completion) / completion : 0.,
			(double) ndp / params.w, hrate, hprobe, hnrate);
		fflush(stdout);
	}

	/*
	 * `red` is the 7-way SUM reduction gathered at the end of the round:
	 *   0: f evaluations   1: collisions   2: probe failures   3: walk robin-hood
	 *   4: walk non-colliding   5: same-value collisions   6: DP failures
	 */
	void end_round(const u64 red[7], const vector<u8> &hll_i)
	{
		double delta = wtime() - round_start;
		u64 N = 1ull << 30;                       /* placeholder, replaced below */
		(void) N;

		ndp_total += ndp;
		ncoll_total += red[1];
		nf_total += red[0];

		for (int k = 0; k < 0x10000; k++)
			if (hll_i[k] > hll[k])
				hll[k] = hll_i[k];

		char hrate[8], hnrate[8];
		human_format((double) red[0] / params.n_walkers / delta, hrate);
		human_format((double) ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);

		printf("\n");
		printf("Round %" PRId64 ".  %.1fs.  #DP %.2f*w (total 2^%.2f).  #coll %.2f*w (total 2^%.2f).  "
		       "Total #f=2^%.3f.  %s #f/s per walker.  node-->%sB/s\n",
			nround, delta,
			(double) ndp / params.w, std::log2((double) ndp_total ? (double) ndp_total : 1.),
			(double) red[1] / params.w, std::log2((double) ncoll_total ? (double) ncoll_total : 1.),
			std::log2((double) nf_total ? (double) nf_total : 1.), hrate, hnrate);

		if (ndp > 0)
			printf("            %.2f%% probe failure.  %.2f%% walk-robinhood.  %.2f%% walk-noncolliding.  "
			       "%.2f%% same-value.  %.2f%% DP failure\n",
				100. * red[2] / ndp, 100. * red[3] / ndp, 100. * red[4] / ndp,
				100. * red[5] / ndp, 100. * red[6] / ndp);

		if (drop_walkerq | drop_out | drop_inserterq | drop_coll)
			printf("            DROPPED  %" PRId64 " walker-queue / %" PRId64 " output-buffer / "
			       "%" PRId64 " inserter-queue / %" PRId64 " collision-queue\n",
				drop_walkerq, drop_out, drop_inserterq, drop_coll);

		u64 E_i = Counters::distinct_collisions_estimation(hll_i);
		u64 E = Counters::distinct_collisions_estimation(hll);
		printf("            #distinct coll (this i / total) %.02f*w / 2^%.2f\n",
			(double) E_i / params.w, std::log2((double) E ? (double) E : 1.));
		printf("\n");
		fflush(stdout);

		nround += 1;
	}

	void done()
	{
		printf("Completed in %.2fs\n", wtime() - start_time);
		fflush(stdout);
	}
};

}
#endif
