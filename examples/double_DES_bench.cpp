#include <mpi.h>

#include "benchmark.hpp"
#include "driver.hpp"
#include "double_DES_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 56;
    u64 seed = 1337;
    std::string engine = "pcs";   // --engine: the search scheme
    mitm::init(argc, argv, opts, ram, n, seed, engine);

    mitm::PRNG prng(seed);
    mitm::DoubleDES_Problem pb(n, prng);
    mitm::benchmark(pb, opts);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
