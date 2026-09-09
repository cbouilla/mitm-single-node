#ifndef MITM_PCS_CONTROL
#define MITM_PCS_CONTROL

#include <mpi.h>

#include "tools.hpp"
#include "pcs/params.hpp"
#include "pcs/shared.hpp"

/*
 * The control channel: a round of PCS ends when enough distinguished points have been found over ALL the
 * nodes, which no node can tell by itself.  So every node's service thread reports its progress to a
 * controller on rank 0, which counts and sends the end of the round back when the count is in -- or as
 * soon as a node says it holds the golden pair, which saves the rest of the round.  Every node then
 * drains and leaves.  The exact numbers are the epilogue's; these are approximate on purpose.
 */

namespace mitm::pcs {

class Control {
	const Params &params;                           /* the round quota, the pacing and the communicator */
	MPI_Request req_end_round = MPI_REQUEST_NULL;   /* every node: the controller's end-of-round */
	MPI_Request req_report = MPI_REQUEST_NULL;      /* rank 0: one progress report as it lands */
	MPI_Request req_solution = MPI_REQUEST_NULL;    /* rank 0: one node's golden pair */
	u64 report_buf[REC_FOUND];                      /* rank 0: where that report lands */
	u64 solution_buf[3];                            /* rank 0: and that pair */
	u64 prev[REC_FOUND];                            /* this node's tallies when it last reported */
	double last_ping = 0;                           /* wtime() when it did */
	bool golden_sent = false;                       /* it has told the controller that it holds a pair */

public:
	u64 reported[REC_FOUND];       /* rank 0: this round's totals as the reports add up.  Approximate */
	u64 total[REC_FOUND];          /* rank 0: the exact all-time sums, from the epilogues */
	RoundStats round;              /* this round's distinct collisions: this node's, then every node's */
	RoundStats all;                /* rank 0: every round's, folded */
	bool round_closed = false;     /* the end-of-round went out: at most one per round */
	bool found = false;            /* a node reported a golden pair: close the round at once */
	double round_start = 0;        /* wtime() when this round began */
	double last_display = 0;       /* wtime() of the last live line */

	Control(const Params &params) : params(params)
	{
		for (int k = 0; k < REC_FOUND; k++) {
			reported[k] = 0;
			total[k] = 0;
			prev[k] = 0;
		}
		MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, params.mpi_comm, &req_end_round);
		if (params.rank != 0)
			return;                /* nobody ever reports to us */
		MPI_Irecv(report_buf, REC_FOUND, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, params.mpi_comm,
		          &req_report);
		MPI_Irecv(solution_buf, 3, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, params.mpi_comm,
		          &req_solution);
	}

	void begin_round()
	{
		for (int k = 0; k < REC_FOUND; k++) {
			reported[k] = 0;
			prev[k] = 0;
		}
		round.reset();
		round_closed = false;
		golden_sent = false;
		round_start = wtime();
		last_ping = round_start;
		last_display = round_start;
	}

	/*
	 * One control turn, from the service thread between two Router turns: nothing here blocks.  `cur` is
	 * this node's record so far -- its tallies then its Router stats -- of which we report the delta.
	 */
	void service(const u64 cur[], Shared &shared)
	{
		int flag = 0;

		MPI_Test(&req_end_round, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			shared.round_over.store(1, std::memory_order_release);
			MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, params.mpi_comm, &req_end_round);
		}

		if (not golden_sent && shared.found.load(std::memory_order_acquire)) {
			MPI_Bsend(shared.golden, 3, MPI_UINT64_T, 0, TAG_SOLUTION, params.mpi_comm);
			golden_sent = true;
		}

		/* report by volume and by time both: on a timer alone a short round overshoots its quota */
		double now = wtime();
		bool due = cur[N_DP] - prev[N_DP] >= params.report_points || now - last_ping >= params.ping_delay;
		if (due && not shared.round_over.load(std::memory_order_relaxed)) {
			u64 msg[REC_FOUND];
			for (int k = 0; k < REC_FOUND; k++) {
				msg[k] = cur[k] - prev[k];
				prev[k] = cur[k];
			}
			last_ping = now;
			MPI_Bsend(msg, REC_FOUND, MPI_UINT64_T, 0, TAG_REPORT, params.mpi_comm);
		}

		if (params.rank != 0)
			return;

		/*
		 * Both receives are drained to empty, which is what keeps a report matched late in our own
		 * drain from being credited to the next round.
		 */
		for (;;) {
			MPI_Test(&req_solution, &flag, MPI_STATUS_IGNORE);
			if (not flag)
				break;
			found = true;
			MPI_Irecv(solution_buf, 3, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_SOLUTION, params.mpi_comm,
			          &req_solution);
		}
		for (;;) {
			MPI_Test(&req_report, &flag, MPI_STATUS_IGNORE);
			if (not flag)
				break;
			for (int k = 0; k < REC_FOUND; k++)
				reported[k] += report_buf[k];
			MPI_Irecv(report_buf, REC_FOUND, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_REPORT, params.mpi_comm,
			          &req_report);
		}

		/* close the round, exactly once: this is where the search advances */
		if (round_closed)
			return;
		if (not found && reported[N_DP] < params.points_per_version)
			return;
		round_closed = true;
		for (int d = 0; d < params.n_nodes; d++)
			MPI_Bsend(NULL, 0, MPI_UINT64_T, d, TAG_END_ROUND, params.mpi_comm);
	}

	/* cancel the posted receives, on the master thread once the team is gone */
	void shutdown()
	{
		if (req_end_round != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_end_round);
			MPI_Wait(&req_end_round, MPI_STATUS_IGNORE);
		}
		if (req_report != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_report);
			MPI_Wait(&req_report, MPI_STATUS_IGNORE);
		}
		if (req_solution != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_solution);
			MPI_Wait(&req_solution, MPI_STATUS_IGNORE);
		}
	}
};

}
#endif
