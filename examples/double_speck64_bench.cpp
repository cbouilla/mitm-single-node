#include <mpi.h>

#include "benchmark.hpp"
#include "driver.hpp"
#include "double_speck64_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 32;
    u64 seed = 1337;
    mitm::init(argc, argv, opts, ram, n, seed);

    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem pb(n, prng);
    mitm::benchmark(pb, opts);
    if (ram > 0) {       // --ram is optional here: it asks for the probe and staging benchmarks too
        mitm::probe_benchmark(pb, opts, ram);
        mitm::staging_benchmark(pb, opts, ram);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
