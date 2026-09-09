#include <mpi.h>

#include "direct.hpp"
#include "driver.hpp"
#include "sha2_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    mitm::Driver drv;
    mitm::init(argc, argv, opts, drv);

    mitm::PRNG prng(drv.seed);
    if (opts.verbose)
        printf("sha2-collision demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, drv.n);

    mitm::SHA2CollisionProblem pb(drv.n, prng);
    mitm::benchmark_and_exit(pb, drv, opts);   // --benchmark: the f/g rate, then exit

    optional<pair<u64, u64>> collision = mitm::direct::collision_search(pb, drv.ram, opts, prng);
    if (opts.verbose) {
        if (collision) {
            auto [x0, x1] = *collision;
            printf("f(%" PRIx64 ") = f(%" PRIx64 ")\n", x0, x1);
        } else {
            printf("Golden collision not found\n");
        }
    }

    MPI_Finalize();
    return collision ? EXIT_SUCCESS : EXIT_FAILURE;
}
