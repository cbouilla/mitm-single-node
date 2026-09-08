#include <mpi.h>

#include "direct.hpp"
#include "driver.hpp"
#include "double_speck64_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 20;          // default problem size (easy)
    u64 seed = 0;        // 0 == draw a fresh one
    std::string engine = "direct";   // --engine: the search scheme
    mitm::init(argc, argv, opts, ram, n, seed, engine);

    mitm::PRNG prng(seed);
    if (opts.verbose)
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n);

    mitm::DoubleSpeck64_Problem Pb(n, prng);
    optional<pair<u64, u64>> claw = mitm::direct::claw_search(Pb, ram, opts, prng);
    if (opts.verbose) {
        if (claw) {
            auto [x0, x1] = *claw;
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        } else {
            printf("Golden claw not found\n");
        }
    }

    MPI_Finalize();
    return claw ? EXIT_SUCCESS : EXIT_FAILURE;
}
