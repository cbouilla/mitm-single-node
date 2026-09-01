#ifndef MITM_PARAMETERS
#define MITM_PARAMETERS

#include <mpi.h>
#include <err.h>
#include <sched.h>
#include <pthread.h>
#include <cmath>
#include <vector>

#include "tools.hpp"
#include "dict.hpp"

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

enum tags {
	TAG_POINTS,        /* bulk distinguished points, node -> node                 */
	TAG_CONTROL,       /* the control channel, both directions (see comm.hpp)     */
};

enum assignment {KEEP_GOING, NEW_VERSION};

enum thread_role {COMM, INSERTER, WALKER};

/* how many u64 per distinguished point on the wire */
static const int DP_WORDS = 3;

/* One distinguished point.  Identical layout on the SPSC queues and on the wire;
   the routing helpers that pick it apart live in comm.hpp. */
struct DP {
	u64 seed;      /* the chain index j it was started from */
	u64 x;         /* the FULL endpoint */
	u64 len;       /* trail length */
};


/******************************* CPU affinity ********************************/

/* the k-th usable CPU of `mask` (k is 0-based); -1 if there are not that many */
static inline int nth_cpu(const cpu_set_t *mask, int k)
{
	for (int cpu = 0; cpu < CPU_SETSIZE; cpu++)
		if (CPU_ISSET(cpu, mask) && (k-- == 0))
			return cpu;
	return -1;
}

/* pin the calling thread to `cpu`.  Returns cpu, or -1 on failure. */
static inline int pin_to_cpu(int cpu)
{
	if (cpu < 0)
		return -1;
	cpu_set_t set;
	CPU_ZERO(&set);
	CPU_SET(cpu, &set);
	if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
		return -1;
	return cpu;
}


/********************************* parameters ********************************/

class Parameters {
public:
	/* MPI topology.  One rank per node, so `rank` IS the node index. */
	MPI_Comm world_comm;
	int rank;
	int n_nodes = 1;                       /* == the number of MPI ranks */

	/* thread layout inside a node (CLI-settable) */
	int walkers_per_node = 0;              /* 0 == fill the inherited affinity mask */
	int inserters_per_node = 1;
	int n_threads;                         /* derived: 1 + walkers_per_node + inserters_per_node */
	int n_walkers;                         /* derived: total walker threads, all nodes */
	int n_inserters;                       /* derived: total dictionary shards, all nodes */

	/* thread pinning */
	bool bind_threads = true;
	cpu_set_t inherited_mask;              /* whatever the launcher gave this rank */
	int n_avail_cpu;

	/* hardware-dependent */
	u64 nbytes_memory = 0;                 /* how much RAM to use on each node */

	/* algorithm parameters */
	double alpha = 2.5;                    /* auto-chosen theta == alpha * sqrt(w/n) */
	double beta = 8;                       /* use each function variant for beta*w DPs */
	double theta = -1;                     /* proportion of distinguished points. -1 == auto */
	u64 multiplier = 0x2545f4914f6cdd1dull;   /* to generate starting points */
	u64 max_versions = 0xffffffffffffffffull; /* how many functions to try before giving up */

	/* SPSC queue capacities, in DPs (rounded up to a power of two internally) */
	size_t walker_queue_capacity = 1024;   /* walker -> comm */
	size_t inserter_queue_capacity = 4096; /* comm -> inserter */

	/* control channel */
	int bsend_slack = 8;                   /* MPI_Bsend slots beyond the n_nodes worst case */

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

	/* other relevant quantities, deduced by finalize() */
	u64 threshold;                         /* any integer less than this is a DP */
	u64 w;                                 /* # slots in the (whole, distributed) dict */
	u64 dp_max_it;                         /* how many iterations to find a DP */
	u64 points_per_version;                /* #DP per version of the function */

	/* utilities */
	bool verbose = 1;                      /* print progress information */

	/*
	 * Discover the topology.  Call this once, right after MPI_Init_thread and before
	 * anything is printed: it is what decides which rank is allowed to talk.
	 */
	void setup(MPI_Comm comm)
	{
		world_comm = comm;
		MPI_Comm_rank(world_comm, &rank);
		MPI_Comm_size(world_comm, &n_nodes);
		verbose = (rank == 0);

		/* what did the launcher actually give us? */
		CPU_ZERO(&inherited_mask);
		if (sched_getaffinity(0, sizeof(inherited_mask), &inherited_mask) != 0)
			err(1, "sched_getaffinity");
		n_avail_cpu = CPU_COUNT(&inherited_mask);

		if (inserters_per_node < 1)
			errx(1, "MPI: at least one inserter thread per node is required");
		if (walkers_per_node == 0)
			walkers_per_node = n_avail_cpu - 1 - inserters_per_node;
		if (walkers_per_node < 1)
			errx(1, "MPI: at least one walker thread per node is required "
			        "(%d CPUs in mask, %d reserved for inserters + comm)",
			        n_avail_cpu, inserters_per_node + 1);

		n_threads = 1 + walkers_per_node + inserters_per_node;
		n_inserters = n_nodes * inserters_per_node;      /* total dictionary shards */
		n_walkers = n_nodes * walkers_per_node;          /* total walker threads */

		if (bind_threads && n_avail_cpu < n_threads)
			warnx("MPI: rank %d has only %d CPUs in its affinity mask for %d threads."
			      "  Did you forget --bind-to none?", rank, n_avail_cpu, n_threads);

		if (verbose) {
			printf("MPI: %d node(s) x (1 comm + %d ins + %d walk) = %d threads/node\n",
				n_nodes, inserters_per_node, walkers_per_node, n_threads);
			printf("MPI: %d dictionary shards, %d walker threads in total\n", n_inserters, n_walkers);
		}
	}

	/* CPU that thread `tid` should run on, or -1 if we cannot place it */
	int cpu_of_thread(int tid) const
	{
		if (!bind_threads || tid >= n_avail_cpu)
			return -1;
		return nth_cpu(&inherited_mask, tid);
	}

	/*
	 * Size the dictionary and choose the difficulty.  Needs the problem dimensions, so
	 * it happens after setup(), once the problem is known.
	 */
	void finalize(int n, int m)
	{
		if (nbytes_memory == 0)
			errx(1, "the amount of RAM to use (per node) must be specified");

		w = PcsDict::get_nslots(nbytes_memory * n_nodes, n_inserters);
		/* auto-choose the difficulty if not set */
		double auto_theta = alpha * std::sqrt((double) w / (1ll << n));
		if (theta < 0) {
			theta = auto_theta;
			if (theta > 1)
				theta = 1;
			if (verbose)
				printf("AUTO-TUNING: setting 1/theta == %.2f\n", 1 / theta);
		} else {
			if (verbose)
				printf("NOTICE: using 1/theta == %.2f vs ``optimal'' 1/theta == %.2f\n", 1/theta, auto_theta);
		}
		threshold = pow(2, m) * theta;
		dp_max_it = 20 / theta;
		points_per_version = beta * w;

		/* display warnings if problematic choices were made */
		if (verbose && theta == 1) {
			printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
			printf("---> zero difficulty (use the naive technique!)\n");
			printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
		}
	}
};


/******************************* benchmarking ********************************/

static void display_stats(u64 N, double start, int vlen, const Parameters &params)
{
	double rate = vlen * N / (wtime() - start);
	double rate_min = rate;
	double rate_max = rate;
	double rate_avg = rate;
	MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, params.world_comm);
	MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, params.world_comm);
	MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
	rate_avg /= params.n_nodes;
	double rate_std = (rate - rate_avg) * (rate - rate_avg);
	MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
	rate_std /= params.n_nodes;
	rate_std = std::sqrt(rate_std);
	if (params.rank == 0) {
		char hmin[8], hmax[8], havg[8], hstd[8];
		human_format(rate_min, hmin);
		human_format(rate_max, hmax);
		human_format(rate_avg, havg);
		human_format(rate_std, hstd);
		printf("Benchmark. f/s: min %s max %s avg %s std %s\n", hmin, hmax, havg, hstd);
	}
}

/* try to iterate for 1s. Return #it/s */
template<typename Problem>
void benchmark(const Problem& pb, const Parameters &params)
{
	if (params.rank == 0)
		printf("Benchmarking scalar implementation (using %d processes)\n", params.n_nodes);

	MPI_Barrier(params.world_comm);

	u64 N = 1ull << 26;
	double start = wtime();
	u64 count = 0;
	for (u64 x = 0; x < N; x++) {
		u64 z = (x & 1) ? pb.f(x) : pb.g(x);
		u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
		int target = ((int) hash) % params.n_inserters;
		if (target == 0)
			count += 1;
	}
	display_stats(N, start, 1, params);

	constexpr int vlen = Problem::vlen;
	if constexpr (vlen > 1) {
		if (params.rank == 0)
			printf("Benchmarking vector implementation (vlen=%d)\n", vlen);

		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 z[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		bool choice[vlen];
		for (int i = 0; i < vlen; i++) {
			choice[i] = i & 1;
			x[i] = i;
		}

		MPI_Barrier(params.world_comm);

		double start = wtime();
		u64 mask = make_mask(pb.n);
		u64 N = 1ull << 20;
		for (u64 i = 0; i < N; i++) {
			pb.vfg(x, choice, z);
			for (int j = 0; j < vlen; j++)
				x[j] = z[j] & mask;
		}
		display_stats(N, start, vlen, params);
	}
}

}
#endif
