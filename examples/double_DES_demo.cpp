#include <mpi.h>

#include "direct.hpp"
#include "driver.hpp"
#include "double_DES_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    mitm::Driver drv;
    mitm::init(argc, argv, opts, drv);

    mitm::PRNG prng(drv.seed);
    if (opts.verbose)
        printf("2DES demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, drv.n);

    mitm::DoubleDES_Problem Pb(drv.n, prng);
    mitm::benchmark_and_exit(Pb, drv, opts);   // --benchmark: the f/g rate, then exit

    optional<pair<u64, u64>> claw = mitm::direct::claw_search(Pb, drv.ram, opts, prng);
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
