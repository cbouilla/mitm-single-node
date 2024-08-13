#ifndef MITM_MPI
#define MITM_MPI

#include <mpi.h>
#include <err.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi_common.hpp"
#include "mpi_sender.hpp"
#include "mpi_receiver.hpp"
#include "mpi_controller.hpp"

namespace mitm {



class MpiEngine {
public:

/* try to iterate for 1s. Return #it/s */
template<typename ConcreteProblem>
static double benchmark(const ConcreteProblem& Pb, MpiParameters &params)
{
	double rate = OpenMPEngine::benchmark(Pb);
	MPI_Allreduce(MPI_IN_PLACE, &rate, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
	return rate;
}

template<typename ConcreteProblem>
static std::tuple<u64,u64,u64> run(const ConcreteProblem& Pb, MpiParameters &params, PRNG &prng)
{
    Dict<std::pair<u64, u64>> dict(params.nbytes_memory);
    params.finalize(Pb.n, dict.n_slots);
    double log2_w = std::log2(dict.n_slots);

    if (params.role == CONTROLLER) {
    	printf("Benchmarking... ");
    	fflush(stdout);
    }
    double it_per_s = benchmark(Pb, params);
    if (params.role == CONTROLLER) {
    	char hitps[8];
    	human_format(it_per_s, hitps);
    	printf("%s iteration/s (using all cores)\n", hitps);
   
    	printf("Starting collision search with seed=%016" PRIx64 ", difficulty=%.2f\n", prng.seed, params.difficulty);
    	printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    	printf("Expected iterations / collision = (2^%0.2f + 2^%.2f) \n", 
    	    Pb.n - params.difficulty - log2_w, 1 + params.difficulty);
    	printf("Expected #iterations = (2^%0.2f + 2^%.2f) \n", 
    	    (Pb.n - 1) + (Pb.n - params.difficulty - log2_w), Pb.n + params.difficulty);
        printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        	params.beta, params.points_per_version, std::log2(params.points_per_version));
    }

    u64 solution[3];
    switch (params.role) {
    case CONTROLLER:
    	solution = controller(Pb, params, prng);
    	break;
    case RECEIVER:
		receiver(Pb, params);
		break;
	case SENDER:
		sender(Pb, params);
	}

	MPI_Bcast(solution, 3, MPI_UINT64_T, 0, params.world_comm);
	return std::tuple(solution[0], solution[1], solution[2]);
}
};


}
#endif