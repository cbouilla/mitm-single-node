#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER
#include <err.h>
#include <vector>

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi_common.hpp"

namespace mitm {



/* Manages send buffers for a collection of receiver processes, with double-buffering */
class SendBuffers : public BaseSendBuffers {
public:
	SendBuffers(const MpiParameters &params) : BaseSendBuffers(params.inter_comm, TAG_POINTS, 3 * params.buffer_capacity) {}
};


template<typename ConcreteProblem>
void sender(const ConcreteProblem& Pb, const MpiParameters &params)
{
	SendBuffers sendbuf(params);
    MpiCounters &ctr = Pb.ctr; 

    u64 steps = 0;
    double last_ping = wtime();

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop

		u64 i = msg[0];
		for (u64 j = msg[1] + 3*params.local_rank;; j += 3*params.n_send) {   // add an odd number to avoid problems mod 2^n...
			steps += 1;

			/* call home? */
            if ((steps % 10000 == 0) && (wtime() - last_ping >= params.ping_delay)) {
				steps = 0;
				last_ping = wtime();
            	u64 stats[4] = {ctr.n_points, ctr.n_dp, (u64) (ctr.send_wait * 1e6), ctr.bytes_sent};
            	MPI_Send(stats, 4, MPI_UINT64_T, 0, TAG_SENDER_CALLHOME, params.world_comm);
            	ctr.reset();

            	int assignment;
            	MPI_Recv(&assignment, 1, MPI_INT, 0, TAG_ASSIGNMENT, params.world_comm, MPI_STATUS_IGNORE);
            	if (assignment == NEW_VERSION) {        /* new broadcast */
            	   	sendbuf.flush(ctr);   
            		break;
            	}
            }

			/* start a new chain from a fresh "random" starting point */
			u64 start = j & Pb.mask;
            auto dp = generate_dist_point(Pb, i, params, start);
            if (not dp) {
                ctr.dp_failure();
                continue;       /* bad chain start */
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);

            u64 hash = (end * 0xdeadbeef) % 0x7fffffff;
            int target_recv = ((int) hash) % params.n_recv;
            sendbuf.push3(start, end, len, target_recv, ctr);
		}
	}
}

}
#endif