#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi/common.hpp"
#include "double_DES_problem.hpp"

u64 seed = 1337;

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    
    mitm::MpiParameters params;
    params.setup(MPI_COMM_WORLD, 0);  // no controller process
    
    mitm::PRNG prng(seed);
    mitm::DoubleDES_Problem pb(56, prng);
    mitm::benchmark(pb, params);

    MPI_Finalize();    
    return EXIT_SUCCESS;
}
