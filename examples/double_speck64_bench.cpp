#include <mpi.h>

#include "benchmark.hpp"
#include "driver.hpp"
#include "double_speck64_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // unused here: the benchmark needs no dictionary
    int n = 32;
    u64 seed = 1337;
    std::string engine = "direct";   // --engine: the search scheme
    mitm::init(argc, argv, opts, ram, n, seed, engine);

    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem pb(n, prng);
    mitm::benchmark(pb, opts);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
