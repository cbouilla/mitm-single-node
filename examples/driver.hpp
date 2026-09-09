#ifndef MITM_DRIVER
#define MITM_DRIVER

#include <mpi.h>
#include <getopt.h>
#include <err.h>
#include <cstdlib>
#include <string>

#include "tools.hpp"
#include "parameters.hpp"
#include "benchmark.hpp"

/*
 * The command line and the MPI startup shared by the drivers.  Part of the examples, not of the
 * library: nothing in ../include depends on it, and a driver of your own may skip it.
 */

namespace mitm {

/*
 * What a driver reads off the command line beside the Options: the defaults below come in, whatever
 * the command line says overwrites them.  Every driver has the same ones, --benchmark included.
 */
struct Driver {
	u64 ram = 0;                       /* --ram: dictionary bytes per node; mandatory for a search */
	int n = 20;                        /* --n: problem size in bits.  Small == easy */
	u64 seed = 0;                      /* --seed: 0 == draw a fresh one from /dev/urandom */
	std::string engine = "direct";     /* --engine: the search scheme */
	bool benchmark = false;            /* --benchmark: measure the problem's f/g rate, then exit */
};

static void usage(const char *argv0)
{
	printf("usage: %s --ram <bytes> [options]\n\n", argv0);
	printf("  --ram B          RAM for the dictionary, per node (accepts 512M, 4G, ...).  MANDATORY\n");
	printf("  --engine S       the search scheme: direct (default).  pcs is being rebuilt on the Router\n");
	printf("  --n BITS         problem size.  Small == easy\n");
	printf("  --seed S         PRNG seed.  0 == draw a fresh one from /dev/urandom\n");
	printf("  --fill F         dictionary fill ratio, entries per round == F * slots.  Default: 0.5\n");
	printf("  --nrounds R      give up after R rounds\n");
	printf("  --benchmark      measure the problem's f/g rate and exit, no search: --ram is then useless\n");
	printf("\n");
	printf("  --producers-per-node W   default: fill the affinity mask\n");
	printf("  --dicts-per-node I       dictionary shards per node, one thread each.  Default: 1\n");
	printf("  --no-bind                do not let the Router pin the threads (needed for 2 ranks per host)\n");
	printf("  --cache-level N          cache level a Router group sits in.  0 == auto.  Default: 0\n");
	printf("  --group N                cores per Router group at most.  Default: 16\n");
	printf("\n");
	printf("  --block N        points per Router block, one block being one message.  Default: 4096\n");
	printf("  --swc N          points a producer accumulates per destination.  0 == auto\n");
	printf("  --n-recv N       MPI receives the Router keeps posted.  Default: 32\n");
	printf("  --inbox N        blocks that may wait for one dict thread.  Default: 64\n");
	printf("  --sweep N        blocks the service thread ships per turn.  Default: 256\n");
	printf("  --credit N       blocks in flight to one peer node at most.  Default: 4\n");
	printf("\n");
	printf("  --difficulty T --alpha A --beta B --dp-len-bits N   (PCS, unused today)\n");
	printf("  --quiet          no progress information\n");
	exit(EXIT_SUCCESS);
}

/* the Router's and the engine's knobs go into `opts`, the driver's own into `drv` */
static void process_command_line_options(int argc, char **argv, Options &opts, Driver &drv)
{
	enum {OPT_ENGINE = 998, OPT_FILL = 999, OPT_BLOCK = 1000, OPT_SWC, OPT_NRECV, OPT_INBOX, OPT_SWEEP,
	      OPT_CREDIT, OPT_GROUP, OPT_CHUNK, OPT_CACHE_LEVEL, OPT_NO_BIND, OPT_DP_LEN_BITS, OPT_BENCHMARK,
	      OPT_QUIET, OPT_HELP};

	struct option longopts[] = {
		{"ram",                required_argument, NULL, 'r'},
		{"engine",             required_argument, NULL, OPT_ENGINE},
		{"n",                  required_argument, NULL, 'n'},
		{"seed",               required_argument, NULL, 's'},
		{"difficulty",         required_argument, NULL, 'd'},
		{"alpha",              required_argument, NULL, 'a'},
		{"beta",               required_argument, NULL, 'b'},
		{"nrounds",            required_argument, NULL, 'o'},
		{"fill",               required_argument, NULL, OPT_FILL},
		{"dp-len-bits",        required_argument, NULL, OPT_DP_LEN_BITS},
		{"producers-per-node", required_argument, NULL, 'W'},
		{"dicts-per-node",     required_argument, NULL, 'I'},
		{"block",              required_argument, NULL, OPT_BLOCK},
		{"swc",                required_argument, NULL, OPT_SWC},
		{"n-recv",             required_argument, NULL, OPT_NRECV},
		{"inbox",              required_argument, NULL, OPT_INBOX},
		{"sweep",              required_argument, NULL, OPT_SWEEP},
		{"credit",             required_argument, NULL, OPT_CREDIT},
		{"group",              required_argument, NULL, OPT_GROUP},
		{"chunk",              required_argument, NULL, OPT_CHUNK},
		{"cache-level",        required_argument, NULL, OPT_CACHE_LEVEL},
		{"no-bind",            no_argument,       NULL, OPT_NO_BIND},
		{"benchmark",          no_argument,       NULL, OPT_BENCHMARK},
		{"quiet",              no_argument,       NULL, OPT_QUIET},
		{"help",               no_argument,       NULL, OPT_HELP},
		{NULL, 0, NULL, 0}
	};

	for (;;) {
		int ch = getopt_long(argc, argv, "", longopts, NULL);
		switch (ch) {
		case -1:                 return;
		case 'r': drv.ram = human_parse(optarg);                             break;
		case OPT_ENGINE:         drv.engine = optarg;                        break;
		case OPT_FILL:           opts.fill = std::stof(optarg);              break;
		case 'n': drv.n = std::stoi(optarg);                                 break;
		case 's': drv.seed = std::stoull(optarg, 0, 0);                      break;
		case 'd': opts.theta = std::stof(optarg);                            break;
		case 'a': opts.alpha = std::stof(optarg);                            break;
		case 'b': opts.beta = std::stof(optarg);                             break;
		case 'o': opts.max_versions = std::stoull(optarg, 0, 0);             break;
		case 'W': opts.producers_per_node = std::stoi(optarg);               break;
		case 'I': opts.dicts_per_node = std::stoi(optarg);                   break;
		case OPT_BLOCK:          opts.router.block_points = human_parse(optarg);  break;
		case OPT_SWC:            opts.router.swc_linesize = human_parse(optarg);  break;
		case OPT_NRECV:          opts.router.n_recv = std::stoi(optarg);      break;
		case OPT_INBOX:          opts.router.inbox_blocks = std::stoi(optarg); break;
		case OPT_SWEEP:          opts.router.sweep_blocks = std::stoi(optarg); break;
		case OPT_CREDIT:         opts.router.credit = std::stoi(optarg);      break;
		case OPT_GROUP:          opts.router.group_size = std::stoi(optarg);  break;
		case OPT_CHUNK:          opts.chunk_size = std::stoull(optarg);       break;
		case OPT_CACHE_LEVEL:    opts.router.cache_level = std::stoi(optarg); break;
		case OPT_DP_LEN_BITS:    opts.dp_lenbits = std::stoi(optarg);         break;
		case OPT_NO_BIND:        opts.router.pin = false;                     break;
		case OPT_BENCHMARK:      drv.benchmark = true;                        break;
		case OPT_QUIET:          opts.verbose = false;                        break;
		case OPT_HELP:           usage(argv[0]);                              break;
		default:
			errx(1, "Unknown option");
		}
	}
}

/*
 * MPI_Init_thread(FUNNELED), the command line, verbose on rank 0 only, and one seed for every rank
 * (drawn by rank 0 and broadcast, if none was given).  `drv` comes in holding the driver's defaults.
 */
static void init(int argc, char **argv, Options &opts, Driver &drv)
{
	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI: this MPI does not provide MPI_THREAD_FUNNELED");

	process_command_line_options(argc, argv, opts, drv);
	if (drv.engine == "pcs")
		errx(1, "--engine pcs: the PCS scheme is being rebuilt on the Router and is not available");
	if (drv.engine != "direct")
		errx(1, "--engine %s: unknown scheme (direct)", drv.engine.c_str());

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	opts.verbose = opts.verbose && (rank == 0);

	if (drv.seed == 0) {
		drv.seed = PRNG::read_urandom();
		MPI_Bcast(&drv.seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	}
}

/*
 * --benchmark: the problem's f/g rate instead of the search.  Every driver calls this once its
 * problem is built; it returns at once unless --benchmark was given, and does not return if it was.
 */
template <class Problem>
static void benchmark_and_exit(const Problem &pb, const Driver &drv, const Options &opts)
{
	if (not drv.benchmark)
		return;
	benchmark(pb, opts);
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}

}
#endif
