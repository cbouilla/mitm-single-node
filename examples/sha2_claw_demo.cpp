#include <mpi.h>

#include "driver.hpp"
#include "sha2_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    mitm::Driver drv;
    mitm::init(argc, argv, opts, drv);

    mitm::PRNG prng(drv.seed);
    if (opts.verbose)
        printf("sha2-claw demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, drv.n);

    mitm::SHA2ClawProblem pb(drv.n, prng);
    mitm::benchmark_and_exit(pb, drv, opts);   // --benchmark: the f/g rate, then exit

    optional<pair<u64, u64>> claw = mitm::claw_search(pb, drv, opts, prng);
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
