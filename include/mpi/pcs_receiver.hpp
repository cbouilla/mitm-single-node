#ifndef MITM_MPI_RECEIVER
#define MITM_MPI_RECEIVER

#include <vector>
#include <mpi.h>

#include "engine_common.hpp"
#include "mpi/common.hpp"

namespace mitm {

template<typename ConcreteProblem>
void receiver(ConcreteProblem& Pb, const MpiParameters &params)
{
    PcsDict dict(params.nbytes_memory / params.recv_per_node);

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop	

		RecvBuffers recvbuf(params.inter_comm, TAG_POINTS, 3 * params.buffer_capacity);
		u64 i = msg[0];
		Pb.n_eval = 0;
		BaseCounters ctr;
		// receive and process data from senders
		for (;;) {
			if (recvbuf.complete())
				break;                      // all senders are done
			auto ready = recvbuf.wait();
			// process incoming buffers of distinguished points
			for (auto it = ready.begin(); it != ready.end(); it++) {
				auto & buffer = **it;
				for (size_t k = 0; k < buffer.size(); k += 3) {
					u64 start = buffer[k];
					u64 end = buffer[k + 1];
					u64 len = buffer[k + 2];
					auto solution = process_distinguished_point(Pb, ctr, dict, i, start, end, len);
					if (solution) {          // call home !
						// maybe save it to a file, just in case
						auto [i, x0, x1] = *solution;
						u64 golden[3] = {i, x0, x1};
						MPI_Send(golden, 3, MPI_UINT64_T, 0, TAG_SOLUTION, params.world_comm);
					}
				}
			}
		}

		// now is a good time to collect stats
		//             #f send     #f recv,      n_collisions,     bytes sent
		u64 imin[4] = {ULLONG_MAX, Pb.n_eval, ctr.n_collisions, ULLONG_MAX};
		u64 imax[4] = {0,          Pb.n_eval, ctr.n_collisions, 0};
		u64 iavg[4] = {0,          Pb.n_eval, ctr.n_collisions, 0};
		MPI_Reduce(imin, NULL, 4, MPI_UINT64_T, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(imax, NULL, 4, MPI_UINT64_T, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(iavg, NULL, 4, MPI_UINT64_T, MPI_SUM, 0, params.world_comm);
		//                send wait recv wait
		double dmin[2] = {HUGE_VAL, recvbuf.waiting_time};
		double dmax[2] = {0,        recvbuf.waiting_time};
		double davg[2] = {0,        recvbuf.waiting_time};
		MPI_Reduce(dmin, NULL, 2, MPI_DOUBLE, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(dmax, NULL, 2, MPI_DOUBLE, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(davg, NULL, 2, MPI_DOUBLE, MPI_SUM, 0, params.world_comm);
		dict.flush();
	}
}

}
#endif
