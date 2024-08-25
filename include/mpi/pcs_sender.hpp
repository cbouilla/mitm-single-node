#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER
#include <err.h>
#include <vector>

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi/common.hpp"

namespace mitm {

template<class ProblemWrapper>
void sender(ProblemWrapper& wrapper, const MpiParameters &params)
{
	u64 mask = make_mask(wrapper.m);
    int jbits = std::log2(10 * params.w) + std::log2(1 / params.theta) + 8;
    u64 jmask = make_mask(jbits);
	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop

    	u64 n_dp = 0;    // #DP found since last report
    	wrapper.n_eval = 0;
		SendBuffers sendbuf(params.inter_comm, TAG_POINTS, 3 * params.buffer_capacity);
    	double last_ping = wtime();
		u64 i = msg[0];
		u64 root_seed = msg[1];
		for (u64 j = params.local_rank ;; j += params.n_send) {

			/* call home? */
            if ((n_dp % 10000 == 9999) && (wtime() - last_ping >= params.ping_delay)) {
				last_ping = wtime();
            	MPI_Send(&n_dp, 1, MPI_UINT64_T, 0, TAG_SENDER_CALLHOME, params.world_comm);
				n_dp = 0;

            	int assignment;
            	MPI_Recv(&assignment, 1, MPI_INT, 0, TAG_ASSIGNMENT, params.world_comm, MPI_STATUS_IGNORE);
            	if (assignment == NEW_VERSION) {        /* new broadcast */
            	   	sendbuf.flush();   
            		break;
            	}
            }

			/* 
			 * start a new chain from a fresh "random" starting point. These are chosen to be
			 * 1) "randomly" spread over [0:2^n]
			 * 2) distinct for each senders
			 * 3) not distinguished
			 */
			if ((j & jmask) != j) {
				printf("j=%" PRIx64 " vs jmask=%" PRIx64 ". n_dp=%" PRIx64 " local_rank=%d, n_send=%d\n", 
					j, jmask, n_dp, params.local_rank, params.n_send);
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			u64 start = (root_seed + j * params.multiplier) & mask;
			if (is_distinguished_point(start, params.threshold))    // refuse to start from a DP
                continue;

            auto dp = generate_dist_point(wrapper, i, params, start);
            if (not dp)
                continue;       /* bad chain start */

			n_dp += 1;
            auto [end, len] = *dp;

            int target_recv = (int) (end % params.n_recv);
            sendbuf.push3(j, end / params.n_recv, len, target_recv);
		}

		// now is a good time to collect stats
		//             #f send,   
		u64 iavg[7] = {wrapper.n_eval, 0, 0, 0, 0, 0, 0};
		MPI_Reduce(iavg, NULL, 7, MPI_UINT64_T, MPI_SUM, 0, params.world_comm);
		//                send wait             recv wait
		double dmin[2] = {sendbuf.waiting_time, HUGE_VAL};
		double dmax[2] = {sendbuf.waiting_time, 0};
		double davg[2] = {sendbuf.waiting_time, 0};
		MPI_Reduce(dmin, NULL, 2, MPI_DOUBLE, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(dmax, NULL, 2, MPI_DOUBLE, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(davg, NULL, 2, MPI_DOUBLE, MPI_SUM, 0, params.world_comm);
	}
}

}
#endif