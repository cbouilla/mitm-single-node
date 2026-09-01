#ifndef MITM_DRIVER
#define MITM_DRIVER

#include <mpi.h>
#include <getopt.h>
#include <err.h>
#include <cstdlib>

#include "tools.hpp"
#include "parameters.hpp"

/*
 * Boilerplate shared by every driver in examples/: the command line is the same for
 * all of them (only the problem changes), and so is the MPI startup sequence.
 */

namespace mitm {

static void usage(const char *argv0)
{
	printf("usage: %s --ram <bytes> [options]\n\n", argv0);
	printf("  --ram B          RAM for the dictionary, per node (accepts 512M, 4G, ...).  MANDATORY\n");
	printf("  --n BITS         problem size.  Small == easy\n");
	printf("  --seed S         PRNG seed.  0 == draw a fresh one from /dev/urandom\n");
	printf("  --difficulty T   proportion theta of distinguished points.  Default: auto\n");
	printf("  --alpha A        auto-tuning: theta = A*sqrt(w/n)\n");
	printf("  --beta B         use each version of the function for B*w distinguished points\n");
	printf("  --nrounds R      give up after R versions of the function\n");
	printf("\n");
	printf("  --walkers-per-node W     default: fill the affinity mask\n");
	printf("  --inserters-per-node I   dictionary shards per node.  Default: 1\n");
	printf("  --no-bind                do not pin threads to CPUs\n");
	printf("\n");
	printf("  --walker-queue N     DPs buffered between a walker and the comm thread\n");
	printf("  --inserter-queue N   DPs buffered between the comm thread and an inserter\n");
	printf("  --coll-queue N       collision candidates buffered for the walkers\n");
	printf("  --coll-per-chunk N   candidates a walker retires per chunk.  0 == drain\n");
	printf("  --buffer N           DPs per node-to-node message\n");
	printf("  --in-buffers N       posted MPI_Irecv slots\n");
	printf("  --chunk N            trail steps between queue / phase checks\n");
	exit(EXIT_SUCCESS);
}

/*
 * `n` and `seed` come in holding the driver's own defaults and are overwritten if the
 * command line says so.
 */
static void process_command_line_options(int argc, char **argv, Parameters &params, int &n, u64 &seed)
{
	enum {OPT_WALKER_QUEUE = 1000, OPT_INSERTER_QUEUE, OPT_COLL_QUEUE, OPT_COLL_PER_CHUNK,
	      OPT_BUFFER, OPT_IN_BUFFERS, OPT_CHUNK, OPT_NO_BIND, OPT_HELP};

	struct option longopts[] = {
		{"ram",                required_argument, NULL, 'r'},
		{"n",                  required_argument, NULL, 'n'},
		{"seed",               required_argument, NULL, 's'},
		{"difficulty",         required_argument, NULL, 'd'},
		{"alpha",              required_argument, NULL, 'a'},
		{"beta",               required_argument, NULL, 'b'},
		{"nrounds",            required_argument, NULL, 'o'},
		{"walkers-per-node",   required_argument, NULL, 'W'},
		{"inserters-per-node", required_argument, NULL, 'I'},
		{"walker-queue",       required_argument, NULL, OPT_WALKER_QUEUE},
		{"inserter-queue",     required_argument, NULL, OPT_INSERTER_QUEUE},
		{"coll-queue",         required_argument, NULL, OPT_COLL_QUEUE},
		{"coll-per-chunk",     required_argument, NULL, OPT_COLL_PER_CHUNK},
		{"buffer",             required_argument, NULL, OPT_BUFFER},
		{"in-buffers",         required_argument, NULL, OPT_IN_BUFFERS},
		{"chunk",              required_argument, NULL, OPT_CHUNK},
		{"no-bind",            no_argument,       NULL, OPT_NO_BIND},
		{"help",               no_argument,       NULL, OPT_HELP},
		{NULL, 0, NULL, 0}
	};

	for (;;) {
		int ch = getopt_long(argc, argv, "", longopts, NULL);
		switch (ch) {
		case -1:                   return;
		case 'r': params.nbytes_memory = human_parse(optarg);            break;
		case 'n': n = std::stoi(optarg);                                 break;
		case 's': seed = std::stoull(optarg, 0, 0);                      break;
		case 'd': params.theta = std::stof(optarg);                      break;
		case 'a': params.alpha = std::stof(optarg);                      break;
		case 'b': params.beta = std::stof(optarg);                       break;
		case 'o': params.max_versions = std::stoull(optarg, 0, 0);       break;
		case 'W': params.walkers_per_node = std::stoi(optarg);           break;
		case 'I': params.inserters_per_node = std::stoi(optarg);         break;
		case OPT_WALKER_QUEUE:   params.walker_queue_capacity = std::stoull(optarg);   break;
		case OPT_INSERTER_QUEUE: params.inserter_queue_capacity = std::stoull(optarg); break;
		case OPT_COLL_QUEUE:     params.coll_queue_capacity = std::stoull(optarg);     break;
		case OPT_COLL_PER_CHUNK: params.coll_per_chunk = std::stoull(optarg);          break;
		case OPT_BUFFER:         params.buffer_capacity = std::stoull(optarg);         break;
		case OPT_IN_BUFFERS:     params.n_in_buffers = std::stoi(optarg);              break;
		case OPT_CHUNK:          params.chunk_size = std::stoull(optarg);              break;
		case OPT_NO_BIND:        params.bind_threads = false;                          break;
		case OPT_HELP:           usage(argv[0]);                                       break;
		default:
			errx(1, "Unknown option");
		}
	}
}

/*
 * Everything a driver does before it can build its problem: start MPI with the
 * threading level the engine requires, read the command line, discover the topology,
 * and make sure every rank uses the same seed -- otherwise the ranks would not even
 * be iterating the same function.
 */
static void init(int argc, char **argv, Parameters &params, int &n, u64 &seed)
{
	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI: this MPI does not provide MPI_THREAD_FUNNELED");

	process_command_line_options(argc, argv, params, n, seed);
	params.setup(MPI_COMM_WORLD);

	if (seed == 0) {
		seed = PRNG::read_urandom();
		MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	}
}

}
#endif
