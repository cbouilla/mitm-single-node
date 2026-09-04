#ifndef MITM_CONTROLLER
#define MITM_CONTROLLER

#include <cmath>
#include <algorithm>
#include <vector>
#include <mpi.h>

#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/*
 * Rank 0's round manager, driven by the comm thread's poll loop (nothing here blocks): it owns the
 * report and solution receives, closes the round with one TAG_END_ROUND per node, and prints.  Exists
 * on every rank so that the node code stays uniform; inert off rank 0.  PROTOCOL.md §2.2 and §4.3.
 */
class Controller {
	MPI_Comm comm;
	MPI_Request req_report = MPI_REQUEST_NULL;        /* rank 0 only */
	MPI_Request req_solution = MPI_REQUEST_NULL;      /* rank 0 only */
	u64 report_buf[N_COUNTERS];
	u64 solution_buf[SOL_NWORDS];

public:
	const Parameters &params;

	u64 nround = 0;                                /* rounds closed so far; the one in progress is the next */
	u64 stop = 0;                                  /* no next round: solution found, or max_versions reached */
	optional<tuple<u64,u64,u64>> solution;         /* (i, x0, x1) */

	/* this round, as the progress reports add up so far.  Approximate (§6 of
	   PROTOCOL.md): it drives the live display and decides when the round closes,
	   nothing else -- the round report prints the exact reduction. */
	u64 reported[N_COUNTERS] = {0};
	bool round_closed = false;                     /* TAG_END_ROUND sent to every node */
	double round_start = 0;                        /* wtime() at begin_round() */
	double last_display = 0;                       /* wtime() of the last live line */

	/* all-time: every round's reduction folded in -- counts summed, HyperLogLog
	   registers maxed.  In-class initialisers: the constructor returns early off rank 0. */
	u64 total[N_COUNTERS] = {0};
	u8 hll[HLL_REGISTERS] = {0};
	double start_time;

	Controller(const Parameters &p)
		: comm(p.mpi_comm), params(p)
	{
		start_time = wtime();
		if (params.rank != 0)
			return;                                /* nobody ever writes to us */
		MPI_Irecv(report_buf, N_COUNTERS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, comm, &req_report);
		MPI_Irecv(solution_buf, SOL_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, comm, &req_solution);
	}

	void begin_round()
	{
		for (int k = 0; k < N_COUNTERS; k++)
			reported[k] = 0;
		round_closed = false;
		round_start = wtime();
		last_display = round_start;
	}

	/*
	 * One turn: drain the solution and report receives to empty, then close the round if it is time
	 * (PROTOCOL.md §2.2).  Draining to empty is what keeps a report matched late in our own drain from
	 * being credited to the next round (§6).
	 */
	void service()
	{
		int flag = 0;

		for (;;) {
			MPI_Test(&req_solution, &flag, MPI_STATUS_IGNORE);
			if (!flag)
				break;
			if (not solution)
				solution = optional(tuple(solution_buf[SOL_I], solution_buf[SOL_X0], solution_buf[SOL_X1]));
			stop = 1;
			MPI_Irecv(solution_buf, SOL_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, comm, &req_solution);
		}

		for (;;) {
			MPI_Test(&req_report, &flag, MPI_STATUS_IGNORE);
			if (!flag)
				break;
			for (int k = 0; k < N_COUNTERS; k++)
				reported[k] += report_buf[k];
			display();
			MPI_Irecv(report_buf, N_COUNTERS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, comm, &req_report);
		}

		/* close the round, exactly once; the search advances here (PROTOCOL.md §2.2) */
		if (not round_closed && (stop || reported[N_DP] >= params.points_per_version)) {
			round_closed = true;
			nround += 1;
			if (nround >= params.max_versions)
				stop = 1;
			for (int r = 0; r < params.n_nodes; r++)
				MPI_Bsend(NULL, 0, MPI_UINT64_T, r, TAG_END_ROUND, comm);
		}
	}

	/* cancel the posted receives; a no-op off rank 0, where none were posted */
	void shutdown()
	{
		if (req_report != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_report);
			MPI_Wait(&req_report, MPI_STATUS_IGNORE);
		}
		if (req_solution != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_solution);
			MPI_Wait(&req_solution, MPI_STATUS_IGNORE);
		}
	}

	/* the live one-line display, from the reports so far; at most twice a second */
	void display()
	{
		double now = wtime();
		if (now - last_display <= 0.5)
			return;
		last_display = now;
		double delta = now - round_start;
		u64 ndp = reported[N_DP];
		double dp_rate = ndp / delta;
		double completion = (double) ndp / params.points_per_version;
		char hrate[8], hnrate[8], hprobe[8];
		human_format(dp_rate / params.theta / params.n_walkers, hrate);
		human_format(ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);
		human_format((double) reported[N_PROBE] / params.n_inserters / delta, hprobe);
		printf("\rRound %" PRId64 ":  %.1fs (%.1f%%, ETA %.1fs).  %.2f*w #DP.  %s #f/s per walker.  "
		       "%s probe/s per inserter.  node-->%sB/s   ",
			nround + 1, delta, 100. * completion,
			(completion > 0) ? delta * (1 - completion) / completion : 0.,
			(double) ndp / params.w, hrate, hprobe, hnrate);
		fflush(stdout);
	}

	/*
	 * The round report: `r` is the exact SUM of the counters over every thread of every node, `hll_round`
	 * the MAX of the HyperLogLog registers; both are folded into the all-time totals.  Printing only:
	 * the round was closed, counted and the next one decided on, in service().
	 */
	void end_round(const u64 r[N_COUNTERS], const u8 hll_round[HLL_REGISTERS])
	{
		double delta = wtime() - round_start;
		for (int k = 0; k < N_COUNTERS; k++)
			total[k] += r[k];
		for (int k = 0; k < HLL_REGISTERS; k++)
			if (hll[k] < hll_round[k])
				hll[k] = hll_round[k];

		u64 ndp = r[N_DP];
		char hrate[8], hnrate[8];
		human_format((double) r[N_EVAL] / params.n_walkers / delta, hrate);
		human_format((double) ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);

		printf("\n");
		printf("Round %" PRId64 ".  %.1fs.  #DP %.2f*w (total 2^%.2f).  #coll %.2f*w (total 2^%.2f).  "
		       "Total #f=2^%.3f.  %s #f/s per walker.  node-->%sB/s\n",
			nround, delta,
			(double) ndp / params.w, std::log2((double) total[N_DP] ? (double) total[N_DP] : 1.),
			(double) r[N_COLLISIONS] / params.w,
			std::log2((double) total[N_COLLISIONS] ? (double) total[N_COLLISIONS] : 1.),
			std::log2((double) total[N_EVAL] ? (double) total[N_EVAL] : 1.), hrate, hnrate);

		if (ndp > 0) {
			double avglen = (double) r[N_POINTS_TRAILS] / ndp;
			printf("            %.2f avg trail length", avglen);
			if (r[N_COLLISIONS] > 0 && avglen > 0)
				printf(" (x%.2f & x%.2f colliding)",
					(double) r[COLLIDING_LEN_MIN] / r[N_COLLISIONS] / avglen,
					(double) r[COLLIDING_LEN_MAX] / r[N_COLLISIONS] / avglen);
			/* BAD_PROBE is an inserter-side counter: its denominator is the probes retired, not the
			   DPs found.  The other four are walker-side, where ndp is right. */
			printf(".  %.2f%% probe failure.  %.2f%% walk-robinhood.  %.2f%% walk-noncolliding.  "
			       "%.2f%% same-value.  %.2f%% DP failure\n",
				r[N_PROBE] ? 100. * r[BAD_PROBE] / r[N_PROBE] : 0., 100. * r[BAD_WALK_ROBINHOOD] / ndp,
				100. * r[BAD_WALK_NONCOLLIDING] / ndp, 100. * r[BAD_COLLISION] / ndp,
				100. * r[BAD_DP] / ndp);

			/*
			 * What the round actually routed.  A round closes on points *found* (service(), against
			 * points_per_version), so a run whose queues overflow burns rounds against a nearly empty
			 * dictionary and still reports them complete: the DPs that reached a shard, and the share of
			 * those found, is the number to read.  PROBLEM.md §2: it is the attack's speed, one for one.
			 */
			char hoff[8], hins[8], hprobe[8];
			human_format((double) ndp / delta, hoff);
			human_format((double) r[N_PROBE] / delta, hins);
			human_format((double) r[N_PROBE] / params.n_inserters / delta, hprobe);
			printf("            ROUTED  %s DP/s found --> %s DP/s inserted (%.2f%% reached a shard,"
			       " %s probe/s per inserter).  dict load %.2f/slot\n",
				hoff, hins, 100. * r[N_PROBE] / ndp, hprobe,
				(double) r[N_PROBE] / params.w);
		}

		if (r[DROP_WALKERQ] | r[DROP_OUT] | r[DROP_INSERTERQ] | r[DROP_COLL])
			printf("            DROPPED  %" PRId64 " walker-queue / %" PRId64 " output-buffer / "
			       "%" PRId64 " inserter-queue / %" PRId64 " collision-queue\n",
				r[DROP_WALKERQ], r[DROP_OUT], r[DROP_INSERTERQ], r[DROP_COLL]);

		u64 E_i = distinct_collisions_estimation(hll_round);
		u64 E = distinct_collisions_estimation(hll);
		printf("            #distinct coll (this i / total) %.02f*w / 2^%.2f\n",
			(double) E_i / params.w, std::log2((double) E ? (double) E : 1.));
		printf("\n");
		fflush(stdout);
	}

	/*
	 * The startup report.  Static: run() prints it before the team, hence the dictionary, exists, so
	 * that the plan is on record even if the allocation fails.
	 */
	static void banner(const Parameters &params, u64 seed)
	{
		char hbuf[8], hdict[8];
		u64 bufbytes = (u64) 2 * params.n_nodes * DP_WORDS * sizeof(u64) * params.buffer_capacity;
		human_format(bufbytes, hbuf);
		human_format(params.w * sizeof(u64), hdict);
		printf("Starting MPI+OpenMP collision search with seed=%016" PRIx64 "\n", seed);
		printf("MPI: %d node(s) x (1 comm + %d ins + %d walk) = %d threads/node\n",
			params.n_nodes, params.inserters_per_node, params.walkers_per_node, params.n_threads);
		printf("MPI: %d dictionary shards, %d walker threads in total\n",
			params.n_inserters, params.n_walkers);
		printf("MPI: rank 0 has %d CPUs in its affinity mask over %d NUMA node(s); threads %s\n",
			params.n_avail_cpu, params.n_numa_nodes, params.bind_threads ? "pinned" : "NOT pinned (--no-bind)");
		int gmin = params.group_size.empty() ? 0 : *std::min_element(params.group_size.begin(),
		                                                             params.group_size.end());
		int gmax = params.group_size.empty() ? 0 : *std::max_element(params.group_size.begin(),
		                                                             params.group_size.end());
		if (params.cache_level_used > 0)
			printf("MPI: %d thread group(s) of %d..%d CPUs over %d L%d domain(s), one inserter each\n",
				params.n_groups, gmin, gmax, params.n_caches, params.cache_level_used);
		else
			printf("MPI: %d thread group(s) of %d..%d CPUs; no cache is shared by several cores,"
			       " so the mask is one domain\n", params.n_groups, gmin, gmax);
		if (params.n_numa_nodes != 1) {
			printf("***** WARNING *****\n");
			printf("---> rank 0 spans %d NUMA nodes: the engine wants ONE MPI RANK PER NUMA NODE\n",
				params.n_numa_nodes);
			printf("---> (mpirun --map-by numa --bind-to numa, or the launcher's equivalent)\n");
			printf("***** WARNING *****\n");
		}
		if (params.bind_threads) {
			printf("MPI: inserters on CPUs");
			for (int tid = 1; tid <= params.inserters_per_node; tid++)
				printf(" %d", params.thread_cpu[tid]);
			printf("\n");
		}
		printf("RAM per node == %sB buffers + dict.  Total dict == %sB (2^%.2f slots)\n",
			hbuf, hdict, std::log2((double) params.w));
		printf("Generating %.1f*w = %" PRIu64 " = 2^%0.2f distinguished points / version\n",
			params.beta, params.points_per_version, std::log2((double) params.points_per_version));
		if (params.theta_auto)
			printf("AUTO-TUNING: setting 1/theta == %.2f\n", 1 / params.theta);
		else
			printf("NOTICE: using 1/theta == %.2f vs ``optimal'' 1/theta == %.2f\n",
				1 / params.theta, 1 / params.auto_theta);
		if (params.theta == 1) {
			printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
			printf("---> zero difficulty (use the naive technique!)\n");
			printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
		}
		fflush(stdout);
	}

	/*
	 * The measured layout, one line per NUMA node: the threads that landed there.  Then the plan, one
	 * line per cache domain, while the groups are few enough to be worth listing.  Static like banner()
	 */
	static void placement(const Parameters &params, const SharedContext &shared)
	{
		std::vector<int> numa_ids;
		for (int tid = 0; tid < params.n_threads; tid++) {
			int id = shared.ctx[tid]->numa_node;
			if (std::find(numa_ids.begin(), numa_ids.end(), id) == numa_ids.end())
				numa_ids.push_back(id);
		}
		std::sort(numa_ids.begin(), numa_ids.end());
		for (size_t k = 0; k < numa_ids.size(); k++) {
			int n_comm = 0, n_ins = 0, n_walk = 0;
			for (int tid = 0; tid < params.n_threads; tid++) {
				if (shared.ctx[tid]->numa_node != numa_ids[k])
					continue;
				if (shared.ctx[tid]->role == COMM)
					n_comm++;
				else if (shared.ctx[tid]->role == INSERTER)
					n_ins++;
				else
					n_walk++;
			}
			printf("NUMA node %d: %s%d ins + %d walk\n", numa_ids[k], n_comm ? "comm + " : "", n_ins, n_walk);
		}
		if (params.cache_level_used > 0 && params.n_groups >= params.n_caches && params.n_groups <= 32)
			for (int c = 0; c < params.n_caches; c++) {
				printf("L%d domain %d: groups", params.cache_level_used, c);
				for (int j = 0; j < params.n_groups; j++)
					if (params.group_cache[j] == c)
						printf(" %d(%d cpu)", j, params.group_size[j]);
				printf("\n");
			}
		else if (params.cache_level_used > 0 && params.n_groups <= 32)
			for (int j = 0; j < params.n_groups; j++)
				printf("group %d: %d cpu, L%d domains %d-%d\n", j, params.group_size[j],
					params.cache_level_used, params.group_cache[j],
					(j + 1 < params.n_groups) ? params.group_cache[j + 1] - 1 : params.n_caches - 1);
		fflush(stdout);
	}

	/* the last line: found, or gave up after nround versions */
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
