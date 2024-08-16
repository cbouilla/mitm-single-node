#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi_common.hpp"
#include "mpi_naive_alltoall.hpp"
#include "mpi_naive_isend.hpp"
#include "double_speck64_problem.hpp"

int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed
bool expensive;

void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    struct option longopts[5] = {
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"recv-per-node", required_argument, NULL, 'e'},
        {"expensive", no_argument, NULL, 'p'},
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
        case 'e':
            params.recv_per_node = std::stoi(optarg);
            break;
        case 'p':
            expensive = 1;
            break;
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mitm::MpiParameters params;
    process_command_line_options(argc, argv, params);
    mitm::PRNG prng(seed);
    if (rank == 0) {
        printf("************************************************************************\n");
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 
    }

    params.setup(MPI_COMM_WORLD, 0);  // no controller process
    mitm::DoubleSpeck64_Problem Pb(n, prng);
    
#if 0
    if (params.verbose) {
        printf("==============================================================\n");
        printf("All-to-all version\n");
        if (expensive) printf("expensive f/g\n");
        printf("==============================================================\n");
    }
    std::vector<std::pair<u64, u64>> claws_alltoall, claws_isend;
    if (expensive)
        claws_alltoall = mitm::naive_mpi_claw_search_alltoall<true>(Pb, params);
    else
        claws_alltoall = mitm::naive_mpi_claw_search_alltoall<false>(Pb, params);
    if (params.verbose)
        for (auto it = claws_alltoall.begin(); it != claws_alltoall.end(); it++) {
            auto [x0, x1] = *it;
            assert(Pb.f(x0) == Pb.g(x1));
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        }
#endif
    if (params.verbose) {
        printf("==============================================================\n");
        printf("Isend version.\n");
        if (expensive) printf("expensive f/g\n");
        printf("==============================================================\n");
    }
    std::vector<std::pair<u64, u64>> claws_isend;
    if (expensive)
        claws_isend = mitm::naive_mpi_claw_search_isend<true>(Pb, params);
    else
        claws_isend = mitm::naive_mpi_claw_search_isend<false>(Pb, params);
    if (params.verbose)
        for (auto it = claws_isend.begin(); it != claws_isend.end(); it++) {
            auto [x0, x1] = *it;
            assert(Pb.f(x0) == Pb.g(x1));
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        }
    assert(claws_isend.size() == 1);
    MPI_Finalize();    
    return EXIT_SUCCESS;
}
