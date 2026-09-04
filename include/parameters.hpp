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
	int cache_level = 0;                      /* cache level a thread group sits in. 0 == auto-detect */
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

	/* collision queues, one per inserter */
	size_t coll_queue_capacity = 8192;     /* candidates buffered by each inserter for its group */
	size_t coll_per_chunk = 0;             /* candidates a walker retires per chunk; 0 == until its
	                                          batch is not full and the queue is empty */

	/* pacing */
	double ping_delay = 0.1;               /* max seconds between reports to the controller */
	int reports_per_round = 4;             /* also report every points_per_version/this DPs,
	                                          so a short round cannot overshoot beta*w */
	size_t chunk_size = 64;                /* vmixf iterations between queue / phase checks */

};


/********************************* thread groups *****************************/

/*
 * Groups are cut out of the CORES of a cache domain, never out of its CPUs: a core's SMT siblings must
 * land in the same group, or the sibling of an inserter would be a thread of somebody else's group.
 *
 * At least as many groups as cache domains: each group sits in one domain, and a domain hosts a share
 * of the groups proportional to its cores, so every group comes out the same size to within one core.
 */
static void split_caches(const std::vector<std::vector<int>> &cache_cores, int n_groups,
                         std::vector<std::vector<int>> &group_cores, std::vector<int> &group_cache)
{
	int n_caches = cache_cores.size();
	std::vector<int> g(n_caches, 1);            /* groups hosted by each domain */
	for (int left = n_groups - n_caches; left > 0; left--) {
		int best = 0;
		double worst = -1;
		for (int c = 0; c < n_caches; c++) {
			double share = (double) cache_cores[c].size() / g[c];
			if (share > worst) {
				worst = share;
				best = c;
			}
		}
		g[best] += 1;                       /* the next group goes where they are largest */
	}

	for (int c = 0; c < n_caches; c++) {
		size_t done = 0;
		for (int k = 0; k < g[c]; k++) {
			size_t left = cache_cores[c].size() - done;
			size_t take = (left + g[c] - k - 1) / (g[c] - k);   /* round up: group 0 also hosts the comm thread */
			group_cores.push_back(std::vector<int>(cache_cores[c].begin() + done,
			                                       cache_cores[c].begin() + done + take));
			group_cache.push_back(c);
			done += take;
		}
	}
}

/*
 * Fewer groups than cache domains: a group must still hold exactly one inserter, so it covers a run of
 * whole domains instead of sitting inside one.  Locality is what is given up; rank 0 warns.
 */
static void merge_caches(const std::vector<std::vector<int>> &cache_cores, int n_groups,
                         std::vector<std::vector<int>> &group_cores, std::vector<int> &group_cache)
{
	int n_caches = cache_cores.size();
	int done = 0;
	for (int j = 0; j < n_groups; j++) {
		int take = (n_caches - done + n_groups - j - 1) / (n_groups - j);   /* round up, as split_caches does */
		std::vector<int> cores;
		for (int c = done; c < done + take; c++)
			cores.insert(cores.end(), cache_cores[c].begin(), cache_cores[c].end());
		group_cores.push_back(cores);
		group_cache.push_back(done);
		done += take;
	}
}

/*
 * The group's emptiest core that still has a free CPU, or -1.  Handing every thread the emptiest core
 * is what puts the comm thread and the inserters on cores of their own, and then fills the siblings
 * they left with walkers: a walker is compute-bound and fills the issue slots a thread waiting on
 * memory leaves idle, whereas two waiting threads on one core would only slow each other down.
 */
static int emptiest_core(const std::vector<int> &cores, const std::vector<std::vector<int>> &core_cpus,
                         const std::vector<size_t> &next_cpu)
{
	int best = -1;
	for (size_t t = 0; t < cores.size(); t++) {
		int k = cores[t];
		if (next_cpu[k] == core_cpus[k].size())
			continue;
		if (best < 0 || next_cpu[k] < next_cpu[best])
			best = k;
	}
	return best;
}


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

	/* thread groups: one per inserter, each inside one cache domain (PROTOCOL.md §1) */
	int cache_level_used;                  /* the level the groups sit in; 0 == none was found */
	int n_caches;                          /* cache domains with at least one CPU in the mask */
	int n_groups;                          /* == inserters_per_node */
	std::vector<int> thread_group;         /* group of each thread, in tid order; a walker's IS its inserter */
	std::vector<int> group_size;           /* CPUs in each group */
	std::vector<int> group_cache;          /* cache domain each group sits in (the first, when it spans) */

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
	 * n, m are the wrapped problem's.  Placement (PROTOCOL.md §1): the mask's CORES are cut into one
	 * group per inserter, each inside one cache domain and all of the same size to within one core; the
	 * comm thread and inserter j then take the emptiest core of group 0 (resp. j), and the walkers fill
	 * the groups evenly, emptiest core first -- so every service thread gets a core of its own and the
	 * SMT siblings it leaves go to walkers.  A walker's group is the inserter it resolves for.
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
		if (bind_threads && inserters_per_node > n_avail_cpu)
			errx(1, "MPI: rank %d: %d shards cannot each own a thread group in a mask of %d CPUs"
			        " (use --no-bind, or fewer shards)", rank, inserters_per_node, n_avail_cpu);
		if (walkers_per_node < inserters_per_node)
			errx(1, "MPI: %d walkers cannot serve %d collision queues: every inserter needs a walker"
			        " of its own group to resolve for it (raise --walkers-per-node, or lower"
			        " --inserters-per-node)", walkers_per_node, inserters_per_node);
		n_threads = 1 + walkers_per_node + inserters_per_node;
		n_inserters = n_nodes * inserters_per_node;
		n_walkers = n_nodes * walkers_per_node;
		n_groups = inserters_per_node;
		if (bind_threads && n_avail_cpu < n_threads)
			warnx("MPI: rank %d has only %d CPUs in its affinity mask for %d threads."
			      "  Did you forget --bind-to none?", rank, n_avail_cpu, n_threads);

		std::vector<int> numa_node_of_cpu;
		std::vector<int> cache_of_cpu;
		std::vector<int> core_of_cpu;
		cache_level_used = topology_of_cpus(mask, cache_level, numa_node_of_cpu, cache_of_cpu, core_of_cpu);
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
			if (cache_level_used > 0 && cache_of_cpu[cpu] < 0) {
				if (rank == 0)
					warnx("MPI: CPU %d of the affinity mask is in no L%d cache hwloc knows:"
					      " the whole mask becomes one domain", cpu, cache_level_used);
				cache_level_used = 0;
			}
		}
		std::sort(numa_ids.begin(), numa_ids.end());
		n_numa_nodes = numa_ids.size();

		/* the mask's CPUs by core, and the cores by cache domain; one domain when no cache is shared */
		std::vector<std::vector<int>> core_cpus;    /* the mask's CPUs of each core, dense index */
		std::vector<int> core_cache;                /* cache domain of each core */
		std::vector<int> dense(CPU_SETSIZE, -1);    /* hwloc's core index -> ours */
		n_caches = 0;
		for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
			if (!CPU_ISSET(cpu, &mask))
				continue;
			int c = (cache_level_used > 0) ? cache_of_cpu[cpu] : 0;
			int hw = core_of_cpu[cpu];
			if (hw < 0 || dense[hw] < 0) {      /* a core not seen yet, or a CPU hwloc puts in none */
				if (hw >= 0)
					dense[hw] = core_cpus.size();
				core_cpus.push_back(std::vector<int>());
				core_cache.push_back(c);
			}
			core_cpus[(hw >= 0) ? dense[hw] : (int) core_cpus.size() - 1].push_back(cpu);
			n_caches = std::max(n_caches, c + 1);
		}
		std::vector<std::vector<int>> cache_cores(n_caches);
		for (size_t k = 0; k < core_cpus.size(); k++)
			cache_cores[core_cache[k]].push_back(k);

		std::vector<std::vector<int>> group_cores;
		if (n_groups >= n_caches) {
			split_caches(cache_cores, n_groups, group_cores, group_cache);
		} else {
			merge_caches(cache_cores, n_groups, group_cores, group_cache);
			if (rank == 0)
				warnx("MPI: %d shards over %d L%d domains: a thread group spans several of them"
				      " (raise --inserters-per-node to %d)", inserters_per_node, n_caches,
				      cache_level_used, n_caches);
		}
		group_size.assign(n_groups, 0);
		int cmin = group_cores[0].size();
		int cmax = 0;
		for (int j = 0; j < n_groups; j++) {
			for (size_t t = 0; t < group_cores[j].size(); t++)
				group_size[j] += core_cpus[group_cores[j][t]].size();
			cmin = std::min(cmin, (int) group_cores[j].size());
			cmax = std::max(cmax, (int) group_cores[j].size());
		}
		if (rank == 0 && cmax > cmin + 1)
			warnx("MPI: thread groups of %d to %d cores: the collision queues are unevenly served",
			      cmin, cmax);

		/*
		 * Every thread takes the emptiest core of its group, service threads first: each of them lands
		 * on a core of its own while cores last, and the siblings they leave are the first the walkers
		 * fill.  A core that runs out is simply skipped.
		 */
		thread_cpu.assign(n_threads, -1);
		thread_group.assign(n_threads, 0);
		std::vector<size_t> next_cpu(core_cpus.size(), 0);   /* core k's next free CPU */
		std::vector<int> group_free(n_groups, 0);            /* CPUs left in each group */
		std::vector<int> n_walk(n_groups, 0);                /* walkers assigned to each group so far */
		int shared = 0;                                      /* service threads that got no core to themselves */
		for (int j = 0; j < n_groups; j++)
			group_free[j] = group_size[j];
		if (bind_threads) {
			int k = emptiest_core(group_cores[0], core_cpus, next_cpu);
			if (k >= 0) {
				thread_cpu[0] = core_cpus[k][next_cpu[k]++];
				group_free[0] -= 1;
			}
		}
		for (int j = 0; j < n_groups; j++) {
			thread_group[1 + j] = j;
			if (not bind_threads)
				continue;
			int k = emptiest_core(group_cores[j], core_cpus, next_cpu);
			if (k < 0)
				continue;
			if (next_cpu[k] > 0)
				shared += 1;                /* the group had fewer cores than service threads */
			thread_cpu[1 + j] = core_cpus[k][next_cpu[k]++];
			group_free[j] -= 1;
		}
		if (shared > 0 && rank == 0)
			warnx("MPI: %d inserter(s) share a core with another service thread: too few cores per"
			      " group for one each", shared);

		/* the emptiest group first, so that no collision queue is left without a walker to drain it */
		for (int s = 0; s < walkers_per_node; s++) {
			int tid = 1 + inserters_per_node + s;
			int best = -1;
			for (int j = 0; bind_threads && j < n_groups; j++)
				if (group_free[j] > 0 && (best < 0 || n_walk[j] < n_walk[best]))
					best = j;
			if (best < 0) {                     /* unpinned, or the mask ran out: spread by index */
				thread_group[tid] = s % n_groups;
				n_walk[s % n_groups] += 1;
				continue;
			}
			int k = emptiest_core(group_cores[best], core_cpus, next_cpu);
			thread_group[tid] = best;
			n_walk[best] += 1;
			thread_cpu[tid] = core_cpus[k][next_cpu[k]++];
			group_free[best] -= 1;
		}

		/* a group whose CPUs ran out borrows a walker from the fullest one: every queue keeps a consumer */
		for (int j = 0; j < n_groups; j++) {
			if (n_walk[j] > 0)
				continue;
			int from = 0;
			for (int k = 0; k < n_groups; k++)
				if (n_walk[k] > n_walk[from])
					from = k;
			for (int tid = 1 + inserters_per_node; tid < n_threads; tid++)
				if (thread_group[tid] == from) {
					thread_group[tid] = j;
					n_walk[from] -= 1;
					n_walk[j] += 1;
					break;
				}
			if (rank == 0)
				warnx("MPI: group %d has no CPU left for a walker: one of group %d resolves for it",
				      j, from);
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
