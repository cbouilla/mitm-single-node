#ifndef MITM_MPI
#define MITM_MPI

#include <mpi.h>
#include <err.h>

#include "common.hpp"
#include "engine_common.hpp"

#include "mpi/common.hpp"
#include "mpi/pcs_sender.hpp"
#include "mpi/pcs_receiver.hpp"
#include "mpi/pcs_controller.hpp"

namespace mitm {

class MpiEngine : Engine {
public:
using Counters = MpiCounters;

template<typename ConcreteProblem>
static tuple<u64,u64,u64> run(const ConcreteProblem& Pb, MpiParameters &params, PRNG &prng)
{
    u64 nslots = PcsDict::get_nslots(params.nbytes_memory / params.recv_per_node) * params.n_recv;
    params.finalize(Pb.n, nslots);

    u64 i, x0, x1;
    switch (params.role) {
    case CONTROLLER:
    	std::tie(i, x0, x1) = controller(Pb, params, prng);
    	break;
    case RECEIVER:
		receiver(Pb, params);
		break;
	case SENDER:
		sender(Pb, params);
	}

	MPI_Bcast(&i, 1, MPI_UINT64_T, 0, params.world_comm);
	MPI_Bcast(&x0, 1, MPI_UINT64_T, 0, params.world_comm);
	MPI_Bcast(&x1, 1, MPI_UINT64_T, 0, params.world_comm);
	return tuple(i, x0, x1);
}
};


}
#endif