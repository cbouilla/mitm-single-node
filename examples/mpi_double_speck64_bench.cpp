#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi_common.hpp"
#include "mpi_naive_alltoall.hpp"
#include "mpi_naive_isend.hpp"
#include "double_speck64_problem.hpp"

u64 seed = 1337;


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    
    mitm::MpiParameters params;
    params.setup(MPI_COMM_WORLD, 0);  // no controller process
    
    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem Pb(32, prng);
    mitm::naive_benchmark(Pb, params);

    MPI_Finalize();    
    return EXIT_SUCCESS;
}
