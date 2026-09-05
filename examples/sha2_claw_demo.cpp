#include <mpi.h>

#include "pcs.hpp"
#include "direct.hpp"
#include "driver.hpp"
#include "sha2_problem.hpp"

int main(int argc, char* argv[])
{
    mitm::Options opts;
    u64 ram = 0;         // RAM per node for the dictionary (--ram, mandatory)
    int n = 20;          // default problem size (easy)
    u64 seed = 0;        // 0 == draw a fresh one
    std::string engine = "pcs";   // --engine: the search scheme
    mitm::init(argc, argv, opts, ram, n, seed, engine);

    mitm::PRNG prng(seed);
    if (opts.verbose)
        printf("sha2-claw demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n);

    mitm::SHA2ClawProblem pb(n, prng);
    optional<pair<u64, u64>> claw;
    if (engine == "direct")
        claw = mitm::direct::claw_search(pb, ram, opts, prng);
    else
        claw = mitm::pcs::claw_search(pb, ram, opts, prng);
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
