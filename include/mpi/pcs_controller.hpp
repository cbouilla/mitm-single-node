#ifndef MITM_MPI_CONTROLLER
#define MITM_MPI_CONTROLLER

#include <cmath>
#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi/common.hpp"

namespace mitm {

/* there is ONE controller process (of global rank 0) */

template<typename ProblemWrapper>
tuple<u64,u64,u64> controller(const ProblemWrapper& wrapper, const MpiParameters &params, PRNG &prng)
{
    printf("Starting MPI collision search with seed=%016" PRIx64 "\n", prng.seed);
    
	char hbsize[8], hdsize[8], htdsize[8];
	u64 bsize_node = 4 * 3 * sizeof(u64) * params.buffer_capacity * params.n_send * params.n_recv / params.n_nodes;
	human_format(bsize_node, hbsize);
	human_format(params.nbytes_memory, hdsize);
	human_format(params.n_nodes * params.nbytes_memory, htdsize);
	double log2_w = std::log2(params.w);
	printf("RAM per node == %sB buffer + %sB dict.  Total dict size == %s (2^%.2f slots)\n", hbsize, hdsize, htdsize, log2_w);

	/* this is quite wrong, actually */
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        	params.beta, params.points_per_version, std::log2(params.points_per_version));

    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
	
	u64 stop = 0;
	u64 nround = 0;
	u64 ndp_total = 0;
	u64 ncoll_total = 0;
	u64 nf_total = 0;
	u64 mask = (1ll << wrapper.pb.m) - 1;
	double start = wtime();

	for (;;) {
        u64 i = prng.rand() & mask;             /* index of families of mixing functions */
       	u64 root_seed = prng.rand();

		/* get data from controller */
		u64 msg[3] = {i, root_seed, stop};
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		
		if (stop)
			break;

		int n_active_senders = params.n_send;
		u64 ndp = 0;                      // #DP found for this i by all senders
		double round_start = wtime();
		double last_display = round_start;
		while (n_active_senders > 0) {
			u64 buffer[10];
			MPI_Status status;
			MPI_Recv(buffer, 10, MPI_UINT64_T, MPI_ANY_SOURCE, MPI_ANY_TAG, params.world_comm, &status);
			switch (status.MPI_TAG) {
				case TAG_SENDER_CALLHOME: {
					ndp += buffer[0];
					int assignment = KEEP_GOING;
					if (stop || ndp >= params.points_per_version) {
						assignment = NEW_VERSION;
						n_active_senders -= 1;
					}
					MPI_Send(&assignment, 1, MPI_INT, status.MPI_SOURCE, TAG_ASSIGNMENT, params.world_comm);

					// verbosity
					double now = wtime();
					if (now - last_display > 0.5) {
						last_display = now;
						double delta = now - round_start;
						double dp_rate = ndp / delta;
						double nf_send_rate = dp_rate / params.theta;
						char hsrate[8], hnrate[8];
						human_format(nf_send_rate / params.n_send, hsrate);
						u64 data_round = ndp * 3 * sizeof(u64) / params.n_nodes;
						human_format(data_round / delta, hnrate);
						double completion = (double) ndp / params.w / params.beta;
						printf("\rRound %" PRId64 ":  %.1fs (%.1f%%, ETA: %.1fs).  %.2f*w #DP.  senders: %s #f/s.  Node-->%sB/s        ",
							nround, delta, 100. * completion, delta / completion, (double) ndp / params.w, hsrate, hnrate);
						fflush(stdout);
					}
					break;
				}

				case TAG_SOLUTION:
					assert(i == buffer[0]);
					solution = optional(tuple(buffer[0], buffer[1], buffer[2]));
					stop = 1;
			}
		}

		// now is a good time to collect and display stats */

		//             #f send, #f recv, collisions, probe_failures, robinhoods, non-colliding, bad_collisions
		u64 iavg[7] = {0, 0, 0, 0, 0, 0, 0};
		MPI_Reduce(MPI_IN_PLACE, iavg, 7, MPI_UINT64_T, MPI_SUM, 0, params.world_comm);
		u64 ncoll = iavg[2];
		ndp_total += ndp;
		ncoll_total += ncoll;
		u64 nf_send = iavg[0];
		u64 nf_recv = iavg[1];
		u64 nf_round = nf_send + nf_recv;
		nf_total += nf_round;

		//                # send wait  #recv wait
		double dmin[2] = {HUGE_VAL,    HUGE_VAL};
		double dmax[2] = {0, 0};
		double davg[2] = {0, 0};
		MPI_Reduce(MPI_IN_PLACE, dmin, 2, MPI_DOUBLE, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(MPI_IN_PLACE, dmax, 2, MPI_DOUBLE, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(MPI_IN_PLACE, davg, 2, MPI_DOUBLE, MPI_SUM, 0, params.world_comm);
		davg[0] /= params.n_send;
		davg[1] /= params.n_recv;

		double delta = wtime() - round_start;

		u64 N = 1ull << wrapper.pb.n;
		char hsrate[8], hrrate[8], hnrate[8];
		human_format(nf_send / params.n_send / delta, hsrate);
		human_format(nf_recv / params.n_recv / delta, hrrate);
		u64 data_round = ndp * 3 * sizeof(u64) / params.n_nodes;
		human_format(data_round / delta, hnrate);
		
		printf("\n");
		printf("Round %" PRId64 " (%.2f*n/w).  %.1fs.  #DP (round / total) %.2f*w / %.2f*n.  #coll (round / total) %.2f*w / %.2f*n.  Total #f=2^%.3f.  node-->%sB/s \n",
			nround, (double) nround * params.w / N, delta, (double) ndp / params.w, (double) ndp_total / N, (double) ncoll / params.w, (double) ncoll_total / N, std::log2(nf_total), hnrate);
		printf("Senders.    Wait == %.2fs / %.2fs (%.1f%%) / %.2fs.  #f == 2^%.2f (%.0f%%).  f/s == %s\n",
                dmin[0], davg[0], 100. * davg[0] / delta, dmax[0], std::log2(nf_send), 100. * nf_send / nf_round, hsrate);
		printf("Receivers.  Wait == %.2fs / %.2fs (%.1f%%) / %.2fs.  #f == 2^%.2f (%.0f%%).  f/s == %s\n",
                dmin[1], davg[1], 100. * davg[1] / delta, dmax[1], std::log2(nf_recv), 100. * nf_recv / nf_round, hrrate);
		printf("            %.2f%% probe failure.  %.2f%% walk-robinhhod.  %.2f%% walk-noncolliding.  %.2f%% same-value\n",
                100. * iavg[3] / ndp, 100. * iavg[4] / ndp, 100. * iavg[5] / ndp, 100. * iavg[6] / ndp);
		printf("\n");
		fflush(stdout);

		nround += 1;
	}
	printf("Completed in %.2fs\n", wtime() - start);

	assert(solution);
	return *solution;
}

}
#endif