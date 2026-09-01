#include <mpi.h>

#include "benchmark.hpp"
#include "driver.hpp"
#include "double_aes_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 32;
    u64 seed = 1337;
    mitm::init(argc, argv, opts, ram, n, seed);

    mitm::PRNG prng(seed);
    mitm::DoubleAES_Problem pb(n, prng);
    mitm::benchmark(pb, opts);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
