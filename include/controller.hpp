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
 * The part of the startup banner every scheme prints: the MPI layout, the thread placement and the
 * dictionary.  A scheme's banner() calls it first, then adds its own lines.
 */
static inline void layout_banner(const Parameters &params)
{
	char hbuf[8], hdict[8];
	u64 bufbytes = (u64) 2 * params.n_nodes * POINT_WORDS * sizeof(u64) * params.buffer_capacity;
	human_format(bufbytes, hbuf);
	human_format(params.w * sizeof(u64), hdict);
	printf("MPI: %d node(s) x (1 comm + %d dict + %d prod) = %d threads/node\n",
		params.n_nodes, params.dicts_per_node, params.producers_per_node, params.n_threads);
	printf("MPI: %d dictionary shards, %d producer threads in total\n",
		params.n_dicts, params.n_producers);
	params.place.report();
	printf("RAM per node == %sB buffers + dict.  Total dict == %sB (2^%.2f slots)\n",
		hbuf, hdict, std::log2((double) params.w));
}


/*
 * Rank 0's round manager, driven by the comm thread's poll loop (nothing here blocks): it owns the
 * report and solution receives, closes the round with one TAG_END_ROUND per node when the scheme says
 * the round is complete or a solution is known, and keeps the tallies.  Exists on every rank so that
 * the node code stays uniform; inert off rank 0.  All printing is the scheme's.  PROTOCOL.md §2.2, §4.3.
 */
template <class Scheme>
class Controller {
	using Params = typename Scheme::Params;
	using RoundStats = typename Scheme::RoundStats;
	static constexpr int N = Scheme::N_COUNTERS;

	MPI_Comm comm;
	MPI_Request req_report = MPI_REQUEST_NULL;        /* rank 0 only */
	MPI_Request req_solution = MPI_REQUEST_NULL;      /* rank 0 only */
	u64 report_buf[N];
	u64 solution_buf[SOL_NWORDS];

public:
	const Params &params;

	u64 nround = 0;                                /* rounds started so far, the one in progress included */
	u64 stop = 0;                                  /* no next round: solution found, or max_versions reached */
	optional<tuple<u64,u64,u64>> solution;         /* (i, x0, x1) */

	/* this round, as the progress reports add up so far.  Approximate (§6 of
	   PROTOCOL.md): it drives the live display and decides when the round closes,
	   nothing else -- the round report prints the exact reduction. */
	u64 reported[N] = {0};
	bool round_closed = false;                     /* TAG_END_ROUND sent to every node */
	double round_start = 0;                        /* wtime() at begin_round() */
	double last_display = 0;                       /* wtime() of the last live line */

	/* all-time: every round's reduction folded in.  In-class initialisers: the
	   constructor returns early off rank 0. */
	u64 total[N] = {0};
	RoundStats all;                                /* the scheme's statistics, folded round after round */
	double start_time;

	Controller(const Params &p)
		: comm(p.mpi_comm), params(p)
	{
		start_time = wtime();
		if (params.rank != 0)
			return;                                /* nobody ever writes to us */
		MPI_Irecv(report_buf, N, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, comm, &req_report);
		MPI_Irecv(solution_buf, SOL_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, comm, &req_solution);
	}

	void begin_round()
	{
		nround += 1;
		for (int k = 0; k < N; k++)
			reported[k] = 0;
		round_closed = false;
		round_start = wtime();
		last_display = round_start;
	}

	/* a solution, off a message or off the epilogue's gather (§4.6): the first one wins and ends the search */
	void record(u64 i, u64 x0, u64 x1)
	{
		if (not solution)
			solution = optional(tuple(i, x0, x1));
		stop = 1;
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
			record(solution_buf[SOL_I], solution_buf[SOL_X0], solution_buf[SOL_X1]);
			MPI_Irecv(solution_buf, SOL_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, comm, &req_solution);
		}

		for (;;) {
			MPI_Test(&req_report, &flag, MPI_STATUS_IGNORE);
			if (!flag)
				break;
			for (int k = 0; k < N; k++)
				reported[k] += report_buf[k];
			display();
			MPI_Irecv(report_buf, N, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, comm, &req_report);
		}

		/* close the round, exactly once; the search advances here (PROTOCOL.md §2.2) */
		if (not round_closed && (stop || Scheme::round_complete(params, reported))) {
			round_closed = true;
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
		Scheme::display(params, reported, now - round_start, nround);
	}

	/*
	 * The round report: `r` is the exact SUM of the counters over every thread of every node, `round`
	 * the scheme's own statistics of the round, reduced likewise; both are folded into the all-time
	 * totals.  Printing only: the round was closed and the next one decided on elsewhere.
	 */
	void end_round(const u64 r[], const RoundStats &round)
	{
		double delta = wtime() - round_start;
		for (int k = 0; k < N; k++)
			total[k] += r[k];
		all.fold(round);
		Scheme::round_report(params, r, total, round, all, delta, nround);
	}
};

}
#endif
