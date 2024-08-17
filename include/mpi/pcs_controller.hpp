#ifndef MITM_MPI_CONTROLLER
#define MITM_MPI_CONTROLLER

#include <cmath>
#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi/common.hpp"

namespace mitm {

/* there is ONE controller process (of global rank 0) */

template<typename ConcreteProblem>
tuple<u64,u64,u64> controller(const ConcreteProblem& Pb, const MpiParameters &params, PRNG &prng)
{
    printf("Starting collision search with seed=%016" PRIx64 ", difficulty=%.2f\n", prng.seed, params.difficulty);
    
	char hbsize[8], hdsize[8], htdsize[8];
	u64 bsize_node = 4 * 3 * sizeof(u64) * params.buffer_capacity * params.n_send * params.n_recv / params.n_nodes;
	human_format(bsize_node, hbsize);
	human_format(params.nbytes_memory, hdsize);
	human_format(params.n_nodes * params.nbytes_memory, htdsize);
	double log2_w = std::log2(params.nslots);
	printf("RAM per node == %sB buffer + %sB dict.  Total dict size == %s (2^%.2f slots)\n", hbsize, hdsize, htdsize, log2_w);

	/* this is quite wrong, actually */
    printf("Expected iterations / collision = (2^%0.2f + 2^%.2f) \n", 
    	    Pb.n - params.difficulty - log2_w, 1 + params.difficulty);
    printf("Expected #iterations = (2^%0.2f + 2^%.2f) \n",
    	    (Pb.n - 1) + (Pb.n - params.difficulty - log2_w), Pb.n + params.difficulty);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        	params.beta, params.points_per_version, std::log2(params.points_per_version));

    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
	
	u64 stop = 0;
	u64 nround = 0;
	u64 ndp_total = 0;
	u64 ncoll_total = 0;
	u64 nf_total = 0;
	double start = wtime();

	for (;;) {
        u64 i = prng.rand() & Pb.mask;             /* index of families of mixing functions */
       	u64 root_seed = prng.rand();

		/* get data from controller */
		u64 msg[3] = {i, root_seed, stop};
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		
		if (stop)
			break;

		int n_active_senders = params.n_send;
		u64 ndp = 0;                      // #DP found for this i by all senders
		double round_start = wtime();
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
					break;
				}

				case TAG_SOLUTION:
					assert(i == buffer[0]);
					solution = optional(tuple(buffer[0], buffer[1], buffer[2]));
					stop = 1;
			}
		}

		// now is a good time to collect and display stats */

		//             #f send,    #f recv,    collisions, bytes sent
		u64 imin[4] = {ULLONG_MAX, ULLONG_MAX, ULLONG_MAX, ULLONG_MAX};
		u64 imax[4] = {0, 0, 0, 0};
		u64 iavg[4] = {0, 0, 0, 0};
		MPI_Reduce(MPI_IN_PLACE, imin, 4, MPI_UINT64_T, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(MPI_IN_PLACE, imax, 4, MPI_UINT64_T, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(MPI_IN_PLACE, iavg, 4, MPI_UINT64_T, MPI_SUM, 0, params.world_comm);
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

		// (min / avg (%% total) / max)
		u64 N = 1ull << Pb.n;
		u64 w = params.nslots;
		char hsrate[8], hrrate[8], hnrate[8];
		human_format(nf_send / params.n_send / delta, hsrate);
		human_format(nf_recv / params.n_recv / delta, hrrate);
		u64 data_round = ndp * 3 * sizeof(u64) / params.n_nodes;
		human_format(data_round / delta, hnrate);
		
		printf("Round %" PRId64 " (%.2f*n/w).  %.1fs.  #DP (round / total) %.2f*w / %.2f*n.  #coll (round / total) %.2f*w / %.2f*n.  Total #f=2^%.3f.  node-->%sB/s \n",
			nround, (double) nround * w / N, delta, (double) ndp / w, (double) ndp_total / N, (double) ncoll / w, (double) ncoll_total / N, std::log2(nf_total), hnrate);
		printf("Senders.    Wait == %.2fs / %.2fs (%.1f%%) / %.2fs.  #f == 2^%.2f (%.0f%%).  f/s == %s\n",
                dmin[0], davg[0], 100. * davg[0] / delta, dmax[0], std::log2(nf_send), 100. * nf_send / nf_round, hsrate);
		printf("Receivers.  Wait == %.2fs / %.2fs (%.1f%%) / %.2fs.  #f == 2^%.2f (%.0f%%).  f/s == %s\n",
                dmin[1], davg[1], 100. * davg[1] / delta, dmax[1], std::log2(nf_recv), 100. * nf_recv / nf_round, hrrate);
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