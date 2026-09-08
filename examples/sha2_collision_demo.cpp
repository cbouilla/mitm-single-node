#include <mpi.h>

#include "direct.hpp"
#include "driver.hpp"
#include "sha2_problem.hpp"

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
        printf("sha2-collision demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n);

    mitm::SHA2CollisionProblem pb(n, prng);
    optional<pair<u64, u64>> collision = mitm::direct::collision_search(pb, ram, opts, prng);
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
