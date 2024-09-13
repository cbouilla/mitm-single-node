#include <cassert>
#include <getopt.h>
#include <err.h>

#include "mitm.hpp"
#include "sequential/pcs_engine.hpp"
#include "double_speck64_problem.hpp"

int n = 20;         // default problem size (easy)
u64 seed = 0x1337;  // default fixed seed

mitm::Parameters process_command_line_options(int argc, char **argv)
{
    struct option longopts[9] = {
        {"ram", required_argument, NULL, 'r'},
        {"difficulty", required_argument, NULL, 'd'},
        {"n", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"nrounds", required_argument, NULL, 'o'},
        {"alpha", required_argument, NULL, 'a'},
        {"beta", required_argument, NULL, 'b'},
        {NULL, 0, NULL, 0}
    };

    mitm::Parameters params;

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return params;
        case 'r':
            params.nbytes_memory = mitm::human_parse(optarg);
            break;
        case 'd':
            params.theta = std::stof(optarg);
            break;
        case 'a':
            params.alpha = std::stof(optarg);
            break;
        case 'b':
            params.beta = std::stof(optarg);
            break;
        case 'n':
            n = std::stoi(optarg);
            break;
        case 's':
            seed = std::stoull(optarg, 0);
            break;
        case 'o':
            params.max_versions = std::stoull(optarg, 0);
            break;            
        default:
            errx(1, "Unknown option %s\n", optarg);
        }
    }
}


int main(int argc, char* argv[])
{
        mitm::Parameters params = process_command_line_options(argc, argv);
        mitm::PRNG prng(seed);
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n); 

        mitm::DoubleSpeck64_Problem Pb(n, prng);            
        auto claw = mitm::claw_search<mitm::ScalarSequentialEngine>(Pb, params, prng);
        if (claw) {
            auto [x0, x1] = *claw;
            printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
        } else {
            printf("Golden collision not found\n");
        }
        return EXIT_SUCCESS;
}
