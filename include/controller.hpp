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
 *
 * It owns the inbound side of rank 0's control channel: one always-posted receive for
 * progress reports (TAG_REPORT) and one for solutions (TAG_SOLUTION), from any node.
 * service() tests them, digests what arrived, and closes the round -- one zero-length
 * TAG_END_ROUND to every node, itself included -- once the reported DP count has
 * reached points_per_version or a solution has come in.  That is the only message the
 * controller ever sends.  The object exists on every rank so that the node code stays
 * uniform, but only rank 0 posts the receives; elsewhere it is inert.
 */
class Controller {
	MPI_Comm comm;
	MPI_Request req_report = MPI_REQUEST_NULL;        /* rank 0 only */
	MPI_Request req_solution = MPI_REQUEST_NULL;      /* rank 0 only */
	u64 report_buf[N_COUNTERS];
	u64 solution_buf[SOL_NWORDS];

public:
	const Parameters &params;

	u64 nround = 0;
	u64 stop = 0;
	optional<tuple<u64,u64,u64>> solution;         /* (i, x0, x1) */

	/* this round, as the progress reports add up so far.  Approximate (§6 of
	   PROTOCOL.md): it drives the live display and decides when the round closes,
	   nothing else -- the round report prints the exact reduction. */
	u64 reported[N_COUNTERS] = {0};
	bool round_closed = false;                     /* TAG_END_ROUND sent to every node */
	double round_start = 0, last_display = 0;

	/* all-time: every round's reduction, merged (counts summed, HyperLogLog maxed) */
	Counters total;
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


	/* Digest one progress report: its deltas are added to the round's tallies. */
	void handle_report(const u64 r[N_COUNTERS])
	{
		for (int k = 0; k < N_COUNTERS; k++)
			reported[k] += r[k];
		display();
	}

	/* A node found the golden pair.  The first one wins, and the search ends. */
	void handle_solution(const u64 s[SOL_NWORDS])
	{
		if (not solution)
			solution = optional(tuple(s[SOL_I], s[SOL_X0], s[SOL_X1]));
		stop = 1;
	}

	/*
	 * Tell every node, ourselves included, that the round is over.  Exactly once per
	 * round: the signal is matched by a receive that names source 0 and this tag, so
	 * successive signals are non-overtaking and a node consumes exactly one per round
	 * (it is its only way out of the steady state) -- a second one would be consumed in
	 * the *next* round and end it at once.  Not to be confused with end_round(), the
	 * statistics, which run once everybody has drained.
	 */
	void close_round()
	{
		round_closed = true;
		for (int r = 0; r < params.n_nodes; r++)
			MPI_Bsend(NULL, 0, MPI_UINT64_T, r, TAG_END_ROUND, comm);
	}

	/*
	 * One turn of the controller: take delivery of everything that has reached rank 0
	 * on the control channel, then decide.  Both inbound messages are one-way and their
	 * relative order is irrelevant: the decision is taken once, after both receives
	 * have been drained.  Each receive is drained to empty so that one call absorbs a
	 * backlog -- this is what keeps a report matched late in our own drain from being
	 * credited to the next round.  Called by the comm thread of rank 0 in every phase of
	 * the round, its own drain included, so late reports and a late solution are
	 * digested (round tallies, and `stop` for the next round header).
	 */
	void service()
	{
		int flag = 0;

		for (;;) {
			MPI_Test(&req_solution, &flag, MPI_STATUS_IGNORE);
			if (!flag)
				break;
			handle_solution(solution_buf);
			MPI_Irecv(solution_buf, SOL_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, comm, &req_solution);
		}

		for (;;) {
			MPI_Test(&req_report, &flag, MPI_STATUS_IGNORE);
			if (!flag)
				break;
			handle_report(report_buf);
			MPI_Irecv(report_buf, N_COUNTERS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, comm, &req_report);
		}

		if (not round_closed && (stop || reported[N_DP] >= params.points_per_version))
			close_round();
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
		printf("\rRound %" PRId64 ":  %.1fs (%.1f%%, ETA %.1fs).  %.2f*w #DP.  %s #f/s per walker.  %s probe/s per inserter.  node-->%sB/s   ",
			nround, delta, 100. * completion,
			(completion > 0) ? delta * (1 - completion) / completion : 0.,
			(double) ndp / params.w, hrate, hprobe, hnrate);
		fflush(stdout);
	}

	/*
	 * `r` is the round's Counters, merged over every thread of every node: the SUM
	 * reduction of the counts and the MAX reduction of the HyperLogLog registers.
	 * Exact, unlike `reported`, so everything printed here comes from it; the
	 * running totals live here too.
	 */
	void end_round(const Counters &r)
	{
		double delta = wtime() - round_start;
		total.merge(r);

		u64 ndp = r.c[N_DP];
		char hrate[8], hnrate[8];
		human_format((double) r.c[N_EVAL] / params.n_walkers / delta, hrate);
		human_format((double) ndp * DP_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);

		printf("\n");
		printf("Round %" PRId64 ".  %.1fs.  #DP %.2f*w (total 2^%.2f).  #coll %.2f*w (total 2^%.2f).  "
		       "Total #f=2^%.3f.  %s #f/s per walker.  node-->%sB/s\n",
			nround, delta,
			(double) ndp / params.w, std::log2((double) total.c[N_DP] ? (double) total.c[N_DP] : 1.),
			(double) r.c[N_COLLISIONS] / params.w,
			std::log2((double) total.c[N_COLLISIONS] ? (double) total.c[N_COLLISIONS] : 1.),
			std::log2((double) total.c[N_EVAL] ? (double) total.c[N_EVAL] : 1.), hrate, hnrate);

		if (ndp > 0) {
			double avglen = (double) r.c[N_POINTS_TRAILS] / ndp;
			printf("            %.2f avg trail length", avglen);
			if (r.c[N_COLLISIONS] > 0 && avglen > 0)
				printf(" (x%.2f & x%.2f colliding)",
					(double) r.c[COLLIDING_LEN_MIN] / r.c[N_COLLISIONS] / avglen,
					(double) r.c[COLLIDING_LEN_MAX] / r.c[N_COLLISIONS] / avglen);
			printf(".  %.2f%% probe failure.  %.2f%% walk-robinhood.  %.2f%% walk-noncolliding.  "
			       "%.2f%% same-value.  %.2f%% DP failure\n",
				100. * r.c[BAD_PROBE] / ndp, 100. * r.c[BAD_WALK_ROBINHOOD] / ndp,
				100. * r.c[BAD_WALK_NONCOLLIDING] / ndp, 100. * r.c[BAD_COLLISION] / ndp,
				100. * r.c[BAD_DP] / ndp);
		}

		if (r.c[DROP_WALKERQ] | r.c[DROP_OUT] | r.c[DROP_INSERTERQ] | r.c[DROP_COLL])
			printf("            DROPPED  %" PRId64 " walker-queue / %" PRId64 " output-buffer / "
			       "%" PRId64 " inserter-queue / %" PRId64 " collision-queue\n",
				r.c[DROP_WALKERQ], r.c[DROP_OUT], r.c[DROP_INSERTERQ], r.c[DROP_COLL]);

		u64 E_i = Counters::distinct_collisions_estimation(r.hll);
		u64 E = Counters::distinct_collisions_estimation(total.hll);
		printf("            #distinct coll (this i / total) %.02f*w / 2^%.2f\n",
			(double) E_i / params.w, std::log2((double) E ? (double) E : 1.));
		printf("\n");
		fflush(stdout);

		nround += 1;
		/* give up after max_versions rounds; the engine then reports "not found" */
		if (nround >= params.max_versions)
			stop = 1;
	}

	/*
	 * The startup report, all of it, in one place.  Static: run() prints it
	 * before the team -- hence the dictionary, and the Controller itself -- exists,
	 * so that the plan is on record even if the allocation fails.
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
