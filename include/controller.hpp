#ifndef MITM_CONTROLLER
#define MITM_CONTROLLER

#include <cmath>
#include <mpi.h>

#include "counters.hpp"
#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/*
 * The controller lives on the comm thread of rank 0, which also has a full node to
 * run.  It is therefore a state object driven by that thread's poll loop rather than
 * a loop of its own: nothing here blocks, and nothing probes.
 */
class Controller {
public:
	const Parameters &params;

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

	Controller(const Parameters &p) : params(p), n_active_nodes(p.n_nodes), hll(0x10000)
	{
		start_time = wtime();
	}

	void banner(const PRNG &prng, u64 w_shard)
	{
		char hbuf[8], hdict[8];
		u64 bufbytes = (u64) 2 * params.n_nodes * DP_WORDS * sizeof(u64) * params.buffer_capacity;
		human_format(bufbytes, hbuf);
		human_format(params.n_nodes * w_shard * sizeof(u64), hdict);
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
	 * `st` is the ST_NWORDS-way SUM reduction of every thread's Counters (see
	 * round_stat in comm.hpp), `hll_round` the MAX-reduction of their HyperLogLog
	 * registers.  Both cover this round only; the running totals live here.
	 */
	void end_round(const u64 st[ST_NWORDS], const vector<u8> &hll_round)
	{
		double delta = wtime() - round_start;

		ndp_total += ndp;
		ncoll_total += st[ST_NCOLL];
		nf_total += st[ST_NEVAL];

		for (int k = 0; k < 0x10000; k++)
			if (hll_round[k] > hll[k])
				hll[k] = hll_round[k];

		char hrate[8], hnrate[8];
		human_format((double) st[ST_NEVAL] / params.n_walkers / delta, hrate);
		human_format((double) ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);

		printf("\n");
		printf("Round %" PRId64 ".  %.1fs.  #DP %.2f*w (total 2^%.2f).  #coll %.2f*w (total 2^%.2f).  "
		       "Total #f=2^%.3f.  %s #f/s per walker.  node-->%sB/s\n",
			nround, delta,
			(double) ndp / params.w, std::log2((double) ndp_total ? (double) ndp_total : 1.),
			(double) st[ST_NCOLL] / params.w, std::log2((double) ncoll_total ? (double) ncoll_total : 1.),
			std::log2((double) nf_total ? (double) nf_total : 1.), hrate, hnrate);

		if (ndp > 0) {
			double avglen = (double) st[ST_NPOINTS_TRAILS] / ndp;
			printf("            %.2f avg trail length", avglen);
			if (st[ST_NCOLL] > 0 && avglen > 0)
				printf(" (x%.2f & x%.2f colliding)",
					(double) st[ST_LEN_MIN] / st[ST_NCOLL] / avglen,
					(double) st[ST_LEN_MAX] / st[ST_NCOLL] / avglen);
			printf(".  %.2f%% probe failure.  %.2f%% walk-robinhood.  %.2f%% walk-noncolliding.  "
			       "%.2f%% same-value.  %.2f%% DP failure\n",
				100. * st[ST_BAD_PROBE] / ndp, 100. * st[ST_BAD_ROBINHOOD] / ndp,
				100. * st[ST_BAD_NONCOLLIDING] / ndp, 100. * st[ST_BAD_COLLISION] / ndp,
				100. * st[ST_BAD_DP] / ndp);
		}

		if (drop_walkerq | drop_out | drop_inserterq | drop_coll)
			printf("            DROPPED  %" PRId64 " walker-queue / %" PRId64 " output-buffer / "
			       "%" PRId64 " inserter-queue / %" PRId64 " collision-queue\n",
				drop_walkerq, drop_out, drop_inserterq, drop_coll);

		u64 E_i = Counters::distinct_collisions_estimation(hll_round);
		u64 E = Counters::distinct_collisions_estimation(hll);
		printf("            #distinct coll (this i / total) %.02f*w / 2^%.2f\n",
			(double) E_i / params.w, std::log2((double) E ? (double) E : 1.));
		printf("\n");
		fflush(stdout);

		nround += 1;
		/* give up after max_versions rounds; the engine then reports "not found" */
		if (nround >= params.max_versions)
			stop = 1;
	}

	void done()
	{
		if (solution)
			printf("Completed in %.2fs\n", wtime() - start_time);
		else
			printf("Gave up after %" PRId64 " versions of the function (%.2fs)\n",
				nround, wtime() - start_time);
		fflush(stdout);
	}
};

}
#endif
