#ifndef MITM_PARAMETERS
#define MITM_PARAMETERS

#include <mpi.h>
#include <err.h>
#include <sched.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>

#include "tools.hpp"

namespace mitm {

/*
 * There is exactly one engine, and it is MPI + OpenMP.
 *
 * Topology: ONE rank per node, MPI_THREAD_FUNNELED.  Inside a rank,
 *     thread 0      == comm (and, on rank 0, the controller),
 *     threads 1..R  == inserters, one dictionary shard each,
 *     the rest      == walkers, which walk trails and resolve collisions.
 * A single rank with one walker and one inserter is the degenerate "sequential" case.
 */

/*
 * Message tags.  TAG_POINTS carries bulk DPs and the end-of-round sentinels between
 * any two nodes; the other three are the control channel: the end-of-round signal goes
 * from rank 0 to every node, reports and solutions from a node to rank 0.  One tag per
 * message kind, so a message is told apart by its envelope and never by its length.
 */
enum tags {TAG_POINTS, TAG_END_ROUND, TAG_REPORT, TAG_SOLUTION};
enum thread_role {COMM, INSERTER, WALKER};

static constexpr int DP_WORDS = 3;     /* how many u64 per distinguished point on the wire */

struct DP {        /* One distinguished point*/
	u64 seed;      /* the chain index j it was started from */
	u64 x;         /* the FULL endpoint */
	u64 len;       /* trail length */
};


/********************************** options **********************************/

struct Options {
	/* the most important option */
	int inserters_per_node = 1;               /* dictionary shards per node. Users SHOULD set this themselves */

	/* other important options */
	double theta = -1;                        /* proportion of distinguished points. -1 == auto */
	bool verbose = true;                      /* print progress information (from rank 0 only) */

	/* --- below this line, the defaults should be just fine */

	/* nodes */
	MPI_Comm mpi_comm = MPI_COMM_WORLD;

	/* thread layout inside a node */
	int walkers_per_node = 0;              /* 0 == fill the inherited affinity mask */
	bool bind_threads = true;              /* pin each thread to a CPU of the mask */

	/* algorithm */
	double alpha = 2.5;                       /* auto-chosen theta == alpha * sqrt(w/n) */
	double beta = 8;                          /* use each function variant for beta*w DPs */
	u64 multiplier = 0x2545f4914f6cdd1dull;   /* to generate starting points */
	u64 max_versions = 0xffffffffffffffffull; /* how many functions to try before giving up */

	/* SPSC queue capacities, in DPs (rounded up to a power of two internally) */
	size_t walker_queue_capacity = 1024;   /* walker -> comm */
	size_t inserter_queue_capacity = 4096; /* comm -> inserter */

	/* control channel */
	int bsend_slack = 8;                   /* MPI_Bsend slots beyond the 2*n_nodes bounded traffic
	                                          (end-of-round signals + sentinels): our reports, our solution */

	/* bulk DP buffering */
	size_t buffer_capacity = 1500;         /* DPs per message (double-buffered per node) */
	int n_in_buffers = 8;                  /* posted MPI_Irecv slots (ANY_SOURCE) */

	/* collision queue */
	size_t coll_queue_capacity = 8192;
	size_t coll_per_chunk = 0;             /* candidates a walker retires per chunk; 0 == drain */

	/* pacing */
	double ping_delay = 0.1;               /* max seconds between reports to the controller */
	int reports_per_round = 4;             /* also report every points_per_version/this DPs,
	                                          so a short round cannot overshoot beta*w */
	size_t chunk_size = 64;                /* vmixf iterations between queue / phase checks */

};


/********************************* parameters ********************************/

/*
 * Everything the engine needs to know, derived ONCE from the options, the RAM budget
 * and the size of the (wrapped) problem: MPI topology, thread layout and placement,
 * dictionary size and difficulty.  Built by run_engine() and never modified
 * afterwards.  Data only: nothing here runs outside the constructor.
 *
 * The user's Options are copied in, so the *resolved* values of walkers_per_node,
 * theta and verbose live here and the caller's object is left alone.
 */
struct Parameters : Options {
	/* MPI topology.  One rank per node, so `rank` IS the node index. */
	int rank;
	int n_nodes;                           /* == the number of MPI ranks */

	/* thread layout */
	int n_threads;                         /* 1 + walkers_per_node + inserters_per_node */
	int n_walkers;                         /* total walker threads, all nodes */
	int n_inserters;                       /* total dictionary shards, all nodes */
	std::vector<int> thread_cpu;           /* CPU for each thread, in tid order; -1 == do not pin */

	/* dictionary */
	u64 nbytes_memory;                     /* RAM budget per node */
	u64 w;                                 /* # slots in the (whole, distributed) dict */
	u64 w_shard;                           /* # slots per inserter thread */
	int jbits;                             /* bits of chain index stored with each DP */

	/* difficulty */
	bool theta_auto;                       /* theta was chosen from alpha */
	double auto_theta;                     /* what alpha chooses, whether or not it was used */
	u64 threshold;                         /* any integer less than this is a DP */
	u64 dp_max_it;                         /* how many iterations to find a DP */
	u64 points_per_version;                /* #DP per version of the function */

	/* n, m are the wrapped problem's: the engine iterates {0,1}^n --> {0,1}^m */
	Parameters(const Options &o, u64 nbytes_memory, int n, int m)
		: Options(o), nbytes_memory(nbytes_memory)
	{
		/* topology */
		MPI_Comm_rank(mpi_comm, &rank);
		MPI_Comm_size(mpi_comm, &n_nodes);
		verbose = verbose && (rank == 0);

		/* what did the launcher actually give us? */
		cpu_set_t mask;
		CPU_ZERO(&mask);
		if (sched_getaffinity(0, sizeof(mask), &mask) != 0)
			err(1, "sched_getaffinity");
		int n_avail_cpu = CPU_COUNT(&mask);

		/* thread layout */
		if (inserters_per_node < 1)
			errx(1, "MPI: at least one inserter thread per node is required");
		if (walkers_per_node == 0)
			walkers_per_node = n_avail_cpu - 1 - inserters_per_node;
		if (walkers_per_node < 1)
			errx(1, "MPI: at least one walker thread per node is required "
			        "(%d CPUs in mask, %d reserved for inserters + comm)",
			        n_avail_cpu, inserters_per_node + 1);
		n_threads = 1 + walkers_per_node + inserters_per_node;
		n_inserters = n_nodes * inserters_per_node;
		n_walkers = n_nodes * walkers_per_node;
		if (bind_threads && n_avail_cpu < n_threads)
			warnx("MPI: rank %d has only %d CPUs in its affinity mask for %d threads."
			      "  Did you forget --bind-to none?", rank, n_avail_cpu, n_threads);

		/* thread placement: the k-th thread goes to the k-th CPU of the mask */
		thread_cpu.assign(n_threads, -1);
		if (bind_threads) {
			int k = 0;
			for (int cpu = 0; cpu < CPU_SETSIZE && k < n_threads; cpu++)
				if (CPU_ISSET(cpu, &mask))
					thread_cpu[k++] = cpu;
		}

		/* dictionary */
		if (nbytes_memory == 0)
			errx(1, "the RAM budget per node must be given (and nonzero)");
		/* 8-byte slots, and a whole number of them per shard */
		w = (nbytes_memory * n_nodes) / sizeof(u64);
		w = (w / n_inserters) * n_inserters;
		if (w == 0)
			errx(1, "RAM budget too small: %" PRIu64 " bytes/node cannot hold one slot per shard",
			     nbytes_memory);
		assert(w % n_inserters == 0);
		w_shard = w / n_inserters;
		jbits = std::log2(10 * w) + 8;

		/* difficulty */
		auto_theta = alpha * std::sqrt(std::ldexp((double) w, -n));
		theta_auto = (theta < 0);
		if (theta_auto)
			theta = std::min(auto_theta, 1.0);
		threshold = pow(2, m) * theta;
		dp_max_it = 20 / theta;
		points_per_version = beta * w;
	}
};

}
#endif
