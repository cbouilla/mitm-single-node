#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "double_speck64_problem.hpp"
#include "mpi/direct.hpp"


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed

void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"threads-per-host", no_argument, NULL, 't'},
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
    mitm::DoubleSpeck64_Problem Pb(n, prng);
    

    // auto console = spdlog::stdout_color_mt("console");
    // spdlog::set_pattern("[%H:%M:%S] [%^%l%$] [tid=%t] %v");
    // spdlog::set_level(spdlog::level::debug); // Set global log level to debug
    // spdlog::debug("This message should be displayed..");    
    

    if (params.verbose) {
        printf("************************************************************************\n");
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 
    }

    vector<pair<u64, u64>> claws = mitm::mpi_direct_claw_search(Pb, params);

    if (params.verbose)
        for (auto it = claws.begin(); it != claws.end(); it++) {
            auto [x0, x1] = *it;
            assert(Pb.f(x0) == Pb.g(x1));
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        }

    assert(claws.size() == 1);
    MPI_Finalize();    
    return EXIT_SUCCESS;
}
