#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi/common.hpp"
#include "double_speck64_problem.hpp"

u64 seed = 1337;

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    
    mitm::MpiParameters params;
    params.setup(MPI_COMM_WORLD);
    
    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem pb(32, prng);
    mitm::benchmark(pb, params);

    MPI_Finalize();    
    return EXIT_SUCCESS;
}
