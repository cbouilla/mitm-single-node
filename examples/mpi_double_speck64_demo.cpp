#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mitm.hpp"
#include "mpi/pcs_engine.hpp"
#include "double_speck64_problem.hpp"


int n = 20;         // default problem size (easy)
u64 seed = 0;       // default random seed


void process_command_line_options(int argc, char **argv, mitm::MpiParameters &params)
{
    enum {OPT_WALKER_QUEUE = 1000, OPT_INSERTER_QUEUE, OPT_COLL_QUEUE,
          OPT_BUFFER, OPT_IN_BUFFERS, OPT_CHUNK, OPT_NO_BIND,
          OPT_COLL_PER_CHUNK};

    struct option longopts[] = {
        {"ram",            required_argument, NULL, 'r'},
        {"n",              required_argument, NULL, 'n'},
        {"seed",           required_argument, NULL, 's'},
        {"difficulty",     required_argument, NULL, 'd'},
        {"alpha",          required_argument, NULL, 'a'},
        {"beta",           required_argument, NULL, 'b'},
        {"walkers-per-node",  required_argument, NULL, 'S'},
        {"inserters-per-node",  required_argument, NULL, 'R'},
        {"walker-queue",     required_argument, NULL, OPT_WALKER_QUEUE},
        {"inserter-queue",     required_argument, NULL, OPT_INSERTER_QUEUE},
        {"coll-queue",     required_argument, NULL, OPT_COLL_QUEUE},
        {"buffer",         required_argument, NULL, OPT_BUFFER},
        {"in-buffers",   required_argument, NULL, OPT_IN_BUFFERS},
        {"chunk",          required_argument, NULL, OPT_CHUNK},
        {"coll-per-chunk", required_argument, NULL, OPT_COLL_PER_CHUNK},
        {"no-bind",        no_argument,       NULL, OPT_NO_BIND},
        {NULL, 0, NULL, 0}
    };

    for (;;) {
        int ch = getopt_long(argc, argv, "", longopts, NULL);
        switch (ch) {
        case -1:
            return;
        case 'r': params.nbytes_memory = mitm::human_parse(optarg); break;
        case 'n': n = std::stoi(optarg); break;
        case 's': seed = std::stoull(optarg, 0); break;
        case 'd': params.theta = std::stof(optarg); break;
        case 'a': params.alpha = std::stof(optarg); break;
        case 'b': params.beta = std::stof(optarg); break;
        case 'S': params.walkers_per_node = std::stoi(optarg); break;
        case 'R': params.inserters_per_node = std::stoi(optarg); break;
        case OPT_WALKER_QUEUE:    params.walker_queue_capacity = std::stoull(optarg); break;
        case OPT_INSERTER_QUEUE:    params.inserter_queue_capacity = std::stoull(optarg); break;
        case OPT_COLL_QUEUE:    params.coll_queue_capacity = std::stoull(optarg); break;
        case OPT_BUFFER:        params.buffer_capacity = std::stoull(optarg); break;
        case OPT_IN_BUFFERS:  params.n_in_buffers = std::stoi(optarg); break;
        case OPT_CHUNK:         params.chunk_size = std::stoull(optarg); break;
        case OPT_COLL_PER_CHUNK: params.coll_per_chunk = std::stoull(optarg); break;
        case OPT_NO_BIND:       params.bind_threads = false; break;
        default:
            errx(1, "Unknown option");
        }
    }
}


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
        errx(1, "MPI: this MPI does not provide MPI_THREAD_FUNNELED");

    mitm::MpiParameters params;
    process_command_line_options(argc, argv, params);
    params.setup(MPI_COMM_WORLD);

    if (seed == 0) {
        seed = mitm::PRNG::read_urandom();
        MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);    // otherwise not everyone evaluates the same function...
    }

    mitm::PRNG prng(seed);
    if (params.rank == 0)
        printf("double-speck64 demo! seed=%016" PRIx64 ", n=%d\n", prng.seed, n);
    mitm::DoubleSpeck64_Problem Pb(n, prng);
    auto claw = mitm::claw_search<mitm::MpiEngine>(Pb, params, prng);
    if (claw && params.rank == 0) {
        auto [x0, x1] = *claw;
        printf("f(%" PRIx64 ") = g(%" PRIx64 ")\n", x0, x1);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
