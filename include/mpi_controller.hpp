#ifndef MITM_MPI_CONTROLLER
#define MITM_MPI_CONTROLLER

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi_common.hpp"

namespace mitm {

/* there is ONE controller process (of global rank 0) */

template<typename ConcreteProblem>
std::tuple<u64,u64,u64> controller(const ConcreteProblem& Pb, const MpiParameters &params, PRNG &prng)
{
	Counters &ctr = Pb.ctr;
    std::optional<std::tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
	u64 stop = 0;
	ctr.ready(Pb.n, params.nslots);
	
	for (;;) {
        u64 i = prng.rand() & Pb.mask;             /* index of families of mixing functions */
       	u64 root_seed = prng.rand();

		/* get data from controller */
		u64 msg[3] = {i, root_seed, stop};
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		
		if (stop)
			break;

		int n_active_senders = params.n_send;
		while (n_active_senders > 0) {
			u64 buffer[10];
			MPI_Status status;
			MPI_Recv(buffer, 10, MPI_UINT64_T, MPI_ANY_SOURCE, MPI_ANY_TAG, params.world_comm, &status);
			switch (status.MPI_TAG) {
				case TAG_SENDER_CALLHOME: {
					ctr.n_points += buffer[0];
					ctr.n_points_trails += buffer[0];
					ctr.n_dp += buffer[1];
					ctr.n_dp_i += buffer[1];
					ctr.send_wait += buffer[2] / 1.0e6;
					ctr.bytes_sent += buffer[3];
					ctr.display();
					int assignment = KEEP_GOING;
					if (stop || ctr.n_dp_i >= params.points_per_version) {
						assignment = NEW_VERSION;
						n_active_senders -= 1;
					}
					MPI_Send(&assignment, 1, MPI_INT, status.MPI_SOURCE, TAG_ASSIGNMENT, params.world_comm);
					break;
				}

				case TAG_RECEIVER_CALLHOME:
					ctr.n_points += buffer[0];
					ctr.n_collisions += buffer[1];
					ctr.n_collisions_i += buffer[1];
					ctr.recv_wait += buffer[2] / 1.0e6;
					ctr.display();
					break;

				case TAG_SOLUTION:
					assert(i == buffer[0]);
					solution = std::make_optional(std::tuple(buffer[0], buffer[1], buffer[2]));
					stop = 1;
			}
		}
		ctr.flush_dict();
	}
	assert(solution);
	return *solution;
}

}
#endif
