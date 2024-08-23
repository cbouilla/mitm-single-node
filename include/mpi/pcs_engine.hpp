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

template<class ProblemWrapper>
static tuple<u64,u64,u64> run(ProblemWrapper& wrapper, MpiParameters &params, PRNG &prng)
{
    u64 i, x0, x1;

    /* safety check: all ranks evaluate the same function */
    u64 test[3];
    u64 mask = make_mask(wrapper.m);
    test[0] = prng.rand() & mask;
    test[1] = prng.rand() & mask;
    test[2] = wrapper.mixf(test[0], test[1]);
    MPI_Bcast(test, 3, MPI_UINT64_T, 0, params.world_comm);
    assert(test[2] == wrapper.mixf(test[0], test[1]));

    /* safety check: w is a multiple of n_recv */
    assert((params.w % params.n_recv) == 0);

    switch (params.role) {
    case CONTROLLER:
    	std::tie(i, x0, x1) = controller(wrapper, params, prng);
    	break;
    case RECEIVER:
		receiver(wrapper, params);
		break;
	case SENDER:
		sender(wrapper, params);
	}

	MPI_Bcast(&i, 1, MPI_UINT64_T, 0, params.world_comm);
	MPI_Bcast(&x0, 1, MPI_UINT64_T, 0, params.world_comm);
	MPI_Bcast(&x1, 1, MPI_UINT64_T, 0, params.world_comm);
	return tuple(i, x0, x1);
}
};


}
#endif