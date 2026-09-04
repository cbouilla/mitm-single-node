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

enum tags {TAG_POINTS, TAG_END_ROUND, TAG_REPORT, TAG_SOLUTION};
enum thread_role {COMM, INSERTER, WALKER};

static constexpr int DP_WORDS = 3;     /* how many u64 per distinguished point on the wire */

struct DP {        /* one distinguished point */
	u64 seed;      /* the chain index j it was started from */
	u64 x;         /* the FULL endpoint */
	u64 len;       /* trail length */
};


/********************************** options **********************************/

struct Options {
	int inserters_per_node = 1;               /* dictionary shards per node. Users SHOULD set this themselves */
	double theta = -1;                        /* proportion of distinguished points. -1 == auto */
	bool verbose = true;                      /* print progress information (from rank 0 only) */

	/* --- below this line, the defaults should be just fine */

	/* nodes */
	MPI_Comm mpi_comm = MPI_COMM_WORLD;       /* one rank per node */

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
	size_t coll_queue_capacity = 8192;     /* collision candidates buffered for the walkers */
	size_t coll_per_chunk = 0;             /* candidates a walker retires per chunk; 0 == drain */

	/* pacing */
	double ping_delay = 0.1;               /* max seconds between reports to the controller */
	int reports_per_round = 4;             /* also report every points_per_version/this DPs,
	                                          so a short round cannot overshoot beta*w */
	size_t chunk_size = 64;                /* vmixf iterations between queue / phase checks */

};


/********************************* parameters ********************************/

/*
 * Everything the engine needs, derived once from (Options, RAM budget, n, m) at the top of run() and
 * const from then on: topology, thread layout and placement, dictionary size, difficulty.  Data only;
 * the Options are copied in, so their resolved values live here.
 */
struct Parameters : Options {
	/* MPI topology.  One rank per node, so `rank` IS the node index. */
	int rank;
	int n_nodes;                           /* == the number of MPI ranks */

	/* thread layout */
	int n_threads;                         /* 1 + walkers_per_node + inserters_per_node */
	int n_walkers;                         /* total walker threads, all nodes */
	int n_inserters;                       /* total dictionary shards, all nodes */
	int n_avail_cpu;                       /* CPUs in the rank's inherited affinity mask */
	int n_numa_nodes;                      /* NUMA nodes with at least one CPU in the mask */
	std::vector<int> thread_cpu;           /* CPU for each thread, in tid order; -1 == do not pin */

	/* dictionary */
	u64 nbytes_memory;                     /* RAM budget per node */
	u64 w;                                 /* # slots in the (whole, distributed) dict */
	u64 w_shard;                           /* # slots per inserter thread */
	int jbits;                             /* bits of chain index stored with each DP */

	/* difficulty */
	bool theta_auto;                       /* theta was chosen from alpha */
	double auto_theta;                     /* what alpha chooses, whether or not it was used */
	u64 threshold;                         /* x is a DP iff x <= threshold */
	u64 dp_max_it;                         /* how many iterations to find a DP */
	u64 points_per_version;                /* #DP per version of the function */

	/*
	 * n, m are the wrapped problem's.  Placement (PROTOCOL.md §1): the comm thread on the first CPU of
	 * the first NUMA node, inserter i and walker s on node i (resp. s) mod n_numa_nodes, next free CPU
	 * each, so the shards are evenly spread only when inserters_per_node is a multiple of n_numa_nodes.
	 */
	Parameters(const Options &o, u64 nbytes_memory, int n, int m)
		: Options(o), nbytes_memory(nbytes_memory)
	{
		MPI_Comm_rank(mpi_comm, &rank);
		MPI_Comm_size(mpi_comm, &n_nodes);
		verbose = verbose && (rank == 0);

		cpu_set_t mask;
		CPU_ZERO(&mask);
		if (sched_getaffinity(0, sizeof(mask), &mask) != 0)
			err(1, "sched_getaffinity");
		n_avail_cpu = CPU_COUNT(&mask);

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

		std::vector<int> numa_node_of_cpu;
		numa_node_of_cpus(mask, numa_node_of_cpu);
		std::vector<int> numa_ids;                  /* distinct NUMA nodes of the mask, ascending */
		for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
			if (!CPU_ISSET(cpu, &mask))
				continue;
			int id = numa_node_of_cpu[cpu];
			if (id < 0 && bind_threads)
				errx(1, "MPI: rank %d: CPU %d of the affinity mask is in no NUMA node hwloc knows"
				        " (use --no-bind)", rank, cpu);
			if (id >= 0 && std::find(numa_ids.begin(), numa_ids.end(), id) == numa_ids.end())
				numa_ids.push_back(id);
		}
		std::sort(numa_ids.begin(), numa_ids.end());
		n_numa_nodes = numa_ids.size();
		if (rank == 0 && n_numa_nodes > 1 && inserters_per_node % n_numa_nodes != 0)
			warnx("MPI: %d shards per node over %d NUMA nodes is not a multiple:"
			      " the dictionary is unevenly spread", inserters_per_node, n_numa_nodes);
		if (bind_threads && n_avail_cpu < 1 + inserters_per_node)
			errx(1, "MPI: rank %d: %d CPUs in the affinity mask cannot pin the comm thread and %d shards"
			        " (use --no-bind, or fewer shards)", rank, n_avail_cpu, inserters_per_node);

		std::vector<std::vector<int>> cpus(n_numa_nodes);   /* the mask's CPUs of each NUMA node */
		for (int cpu = 0; cpu < CPU_SETSIZE; cpu++)
			if (CPU_ISSET(cpu, &mask) && numa_node_of_cpu[cpu] >= 0) {
				int k = std::find(numa_ids.begin(), numa_ids.end(), numa_node_of_cpu[cpu]) - numa_ids.begin();
				cpus[k].push_back(cpu);
			}
		std::vector<size_t> next(n_numa_nodes, 0);          /* cpus[k][next[k]]: node k's next free CPU */
		thread_cpu.assign(n_threads, -1);
		if (bind_threads) {
			for (int tid = 0; tid < n_threads; tid++) {
				int k = 0;
				if (tid >= 1 && tid <= inserters_per_node)
					k = (tid - 1) % n_numa_nodes;
				else if (tid > inserters_per_node)
					k = (tid - 1 - inserters_per_node) % n_numa_nodes;
				int tries = 0;
				while (tries < n_numa_nodes && next[k] == cpus[k].size()) {
					k = (k + 1) % n_numa_nodes;
					tries++;
				}
				if (tries == n_numa_nodes)
					break;                              /* every CPU of the mask is taken */
				thread_cpu[tid] = cpus[k][next[k]];
				next[k]++;
			}
		}

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
