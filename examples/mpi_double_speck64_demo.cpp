#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mitm.hpp"
#include "mpi/pcs_engine.hpp"
#include "double_speck64_problem.hpp"


int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed



mitm::Parameters process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"ram", required_argument, NULL, 'r'},
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"recv-per-node", required_argument, NULL, 'e'},
        {NULL, 0, NULL, 0}
    };

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return params;
        case 'r':
            params.nbytes_memory = mitm::human_parse(optarg);
            break;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 's':
            seed = std::stoull(optarg, 0);
            break;
        case 'e':
            params.recv_per_node = std::stoi(optarg);
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    mitm::MpiParameters params;
    process_command_line_options(argc, argv, params);
    params.setup(MPI_COMM_WORLD);

    mitm::PRNG prng(seed);
    if (params.role == mitm::CONTROLLER)
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 
    mitm::DoubleSpeck64_Problem Pb(n, prng);
    auto claw = mitm::claw_search<mitm::MpiEngine>(Pb, params, prng);
    if (params.role == mitm::CONTROLLER)
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", claw.first, claw.second);
    
    MPI_Finalize();    
    return EXIT_SUCCESS;
}
