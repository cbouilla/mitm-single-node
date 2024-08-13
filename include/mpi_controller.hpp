#ifndef MITM_MPI_CONTROLLER
#define MITM_MPI_CONTROLLER

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi_common.hpp"

namespace mitm {

template<typename ConcreteProblem>
std::tuple<u64,u64,u64> controller(const ConcreteProblem& Pb, const MpiParameters &params, PRNG &prng)
{
	Counters &ctr = Pb.ctr;
    std::optional<std::tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
	u64 i = 0;                     /* index of families of mixing functions */
	u64 stop = 0;

	for (;;) {
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
				case TAG_SENDER_CALLHOME:
					// collect stats
					ctr.n_dp_i += buffer[0];
					int assignment = KEEP_GOING;
					if (stop || ctr.n_dp_i >= params.points_per_version) {
						assignment = NEW_VERSION;
						n_active_senders -= 1;
					}
					MPI_Send(&assignment, 1, MPI_INT, status.MPI_SOURCE, TAG_ASSIGNMENT, params.world_comm);
					break;
				case TAG_RECEIVER_CALLHOME:
					// collect stats
					break;
				case TAG_SOLUTION:
					solution = std::make_optional(std::tuple(buffer[0], buffer[1], buffer[2]));
					printf("Got solution, trying to stop everything\n");
					stop = 1;
			}
			ctr.display();
		}
	}
	assert(solution);
	return *solution;
}

}
#endif
