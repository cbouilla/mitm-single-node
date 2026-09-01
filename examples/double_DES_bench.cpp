#include <mpi.h>

#include "parameters.hpp"
#include "driver.hpp"
#include "double_DES_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Parameters params;
    int n = 56;
    u64 seed = 1337;
    mitm::init(argc, argv, params, n, seed);

    mitm::PRNG prng(seed);
    mitm::DoubleDES_Problem pb(n, prng);
    mitm::benchmark(pb, params);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
