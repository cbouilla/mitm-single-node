#include <mpi.h>

#include "mitm.hpp"
#include "driver.hpp"
#include "double_aes_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 20;          // default problem size (easy)
    u64 seed = 0;        // 0 == draw a fresh one
    mitm::init(argc, argv, opts, ram, n, seed);

    mitm::PRNG prng(seed);
    if (opts.verbose)
        printf("double-aes demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n);

    mitm::DoubleAES_Problem Pb(n, prng);
    auto claw = mitm::claw_search(Pb, ram, opts, prng);
    if (opts.verbose) {
        if (claw) {
            auto [x0, x1] = *claw;
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        } else {
            printf("Golden claw not found\n");
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
