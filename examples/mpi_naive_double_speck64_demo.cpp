#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "double_speck64_problem.hpp"
#include "mpi/direct_alltoall.hpp"
#include "mpi/benchmark.hpp"

int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed
bool bench = 0;

void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"threads-per-host", required_argument, NULL, 't'},
        {"benchmark", no_argument, NULL, 'b'},
        {NULL, 0, NULL, 0}
    };

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 's':
            seed = std::stoull(optarg, 0);
            break;
        case 't':
            params.n_threads = std::stoi(optarg);
            break;
        case 'b':
            bench = 1;
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    assert(provided >= MPI_THREAD_FUNNELED);

    mitm::MpiParameters params(MPI_COMM_WORLD);
    process_command_line_options(argc, argv, params);
    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem pb(n, prng);
    
    if (params.verbose) {
        println("************************************************************************");
        println("double-speck64 demo! seed={:016x}, n={}", prng.seed, n); 
    }

    if (bench)
        mitm::benchmark(pb, params);

    vector<pair<u64, u64>> claws = mitm::mpi_direct_claw_search(pb, params);

    if (params.verbose)
        for (auto it = claws.begin(); it != claws.end(); it++) {
            auto [x0, x1] = *it;
            assert(pb.f(x0) == pb.g(x1));
            println("f({:x}) = g({:x})", x0, x1);
        }

    assert(claws.size() == 1);
    MPI_Finalize();    
    return EXIT_SUCCESS;
}
