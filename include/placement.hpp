#ifndef MITM_PLACEMENT
#define MITM_PLACEMENT

#include <err.h>
#include <sched.h>
#include <cstdio>
#include <algorithm>
#include <vector>

#include "router/topology.hpp"

/*
 * Where the threads of a rank go, and nothing else.  On top of the shared hwloc primitives (topology.hpp),
 * the partition of the affinity mask into thread groups, the CPU each thread is pinned to, and the reports
 * of both the plan and what the kernel actually did.  PROTOCOL.md §1 is the specification.
 */

namespace mitm {

enum thread_role {COMM, DICT, PRODUCER};


/********************************* thread groups *****************************/

/*
 * Groups are cut out of the CORES of a cache domain, never out of its CPUs: a core's SMT siblings must
 * land in the same group, or the sibling of a dict thread would be a thread of somebody else's group.
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
 * Fewer groups than cache domains: a group must still hold exactly one dict thread, so it covers a run of
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

/********************************* the placement *****************************/

/*
 * The rank's thread layout, derived once from its affinity mask and const from then on.  Data and its
 * two reports; the engine reads `thread_cpu` and `thread_group` and nothing else.
 */
struct Placement {
	const int rank;                        /* whose affinity mask this describes */
	const int n_groups;                    /* == dicts_per_node: one group per dict thread */
	const bool bind;                       /* pin, or leave every thread where the launcher put it */
	int n_producers;                       /* producers_per_node, resolved: 0 asked to fill the mask */

	/* what the mask offers */
	int n_avail_cpu;                       /* CPUs in it */
	int n_numa_nodes;                      /* NUMA nodes with at least one of them */
	int cache_level;                       /* the level the groups sit in; 0 == none was found */
	int n_caches;                          /* cache domains with at least one of them */

	/* the plan */
	std::vector<int> thread_cpu;           /* CPU for each thread, in tid order; -1 == do not pin */
	std::vector<int> thread_group;         /* group of each thread; a producer's IS its dict thread */
	std::vector<int> group_size;           /* CPUs in each group */
	std::vector<int> group_cache;          /* cache domain of each group (the first, when it spans) */

	/*
	 * The mask's CORES are cut into one group per dict thread, each inside one cache domain and all of the
	 * same size to within one core; the comm thread and dict thread j then take the emptiest core of group
	 * 0 (resp. j), and the producers fill the groups evenly, emptiest core first -- so every service
	 * thread gets a core of its own and the SMT siblings it leaves go to producers.  `producers == 0` asks
	 * for as many as the mask holds.  `producer_per_dict` is the scheme's: PCS needs a producer in every
	 * group to drain its dict thread's collision queue, the direct scheme does not.
	 */
	Placement(int rank, int dicts, int producers, bool bind, int level, bool producer_per_dict)
		: rank(rank), n_groups(dicts), bind(bind), n_producers(producers)
	{
		cpu_set_t mask;
		CPU_ZERO(&mask);
		if (sched_getaffinity(0, sizeof(mask), &mask) != 0)
			err(1, "sched_getaffinity");
		n_avail_cpu = CPU_COUNT(&mask);

		if (n_groups < 1)
			errx(1, "MPI: at least one dict thread per node is required");
		if (n_producers == 0)
			n_producers = n_avail_cpu - 1 - n_groups;
		if (n_producers < 1)
			errx(1, "MPI: at least one producer thread per node is required "
			        "(%d CPUs in mask, %d reserved for dict threads + comm)", n_avail_cpu, n_groups + 1);
		if (bind && n_groups > n_avail_cpu)
			errx(1, "MPI: rank %d: %d shards cannot each own a thread group in a mask of %d CPUs"
			        " (use --no-bind, or fewer shards)", rank, n_groups, n_avail_cpu);
		if (producer_per_dict && n_producers < n_groups)
			errx(1, "MPI: %d producers cannot serve %d collision queues: every dict thread needs a producer"
			        " of its own group to resolve for it (raise --producers-per-node, or lower"
			        " --dicts-per-node)", n_producers, n_groups);
		int n_threads = 1 + n_groups + n_producers;
		if (bind && n_avail_cpu < n_threads)
			warnx("MPI: rank %d has only %d CPUs in its affinity mask for %d threads."
			      "  Did you forget --bind-to none?", rank, n_avail_cpu, n_threads);

		std::vector<int> numa_node_of_cpu;
		std::vector<int> cache_of_cpu;
		std::vector<int> core_of_cpu;
		cache_level = topology_of_cpus(mask, level, numa_node_of_cpu, cache_of_cpu, core_of_cpu);
		std::vector<int> numa_ids;                  /* distinct NUMA nodes of the mask, ascending */
		for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
			if (!CPU_ISSET(cpu, &mask))
				continue;
			int id = numa_node_of_cpu[cpu];
			if (id < 0 && bind)
				errx(1, "MPI: rank %d: CPU %d of the affinity mask is in no NUMA node hwloc knows"
				        " (use --no-bind)", rank, cpu);
			if (id >= 0 && std::find(numa_ids.begin(), numa_ids.end(), id) == numa_ids.end())
				numa_ids.push_back(id);
			if (cache_level > 0 && cache_of_cpu[cpu] < 0) {
				if (rank == 0)
					warnx("MPI: CPU %d of the affinity mask is in no L%d cache hwloc knows:"
					      " the whole mask becomes one domain", cpu, cache_level);
				cache_level = 0;
			}
		}
		std::sort(numa_ids.begin(), numa_ids.end());
		n_numa_nodes = numa_ids.size();

		/* the mask's CPUs by core, and the cores by cache domain; one domain when no cache is shared */
		std::vector<std::vector<int>> core_cpus;    /* the mask's CPUs of each core, dense index */
		std::vector<std::vector<int>> cache_cores;  /* the dense core indices of each cache domain */
		n_caches = cores_by_domain(mask, cache_level, cache_of_cpu, core_of_cpu, core_cpus, cache_cores);

		std::vector<std::vector<int>> group_cores;
		if (n_groups >= n_caches) {
			split_caches(cache_cores, n_groups, group_cores, group_cache);
		} else {
			merge_caches(cache_cores, n_groups, group_cores, group_cache);
			if (rank == 0)
				warnx("MPI: %d shards over %d L%d domains: a thread group spans several of them"
				      " (raise --dicts-per-node to %d)", n_groups, n_caches, cache_level, n_caches);
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
		if (producer_per_dict && rank == 0 && cmax > cmin + 1)
			warnx("MPI: thread groups of %d to %d cores: the collision queues are unevenly served",
			      cmin, cmax);

		/*
		 * Every thread takes the emptiest core of its group, service threads first: each of them lands
		 * on a core of its own while cores last, and the siblings they leave are the first the producers
		 * fill.  A core that runs out is simply skipped.
		 */
		thread_cpu.assign(n_threads, -1);
		thread_group.assign(n_threads, 0);
		std::vector<size_t> next_cpu(core_cpus.size(), 0);   /* core k's next free CPU */
		std::vector<int> group_free(n_groups, 0);            /* CPUs left in each group */
		std::vector<int> n_walk(n_groups, 0);                /* producers assigned to each group so far */
		int shared = 0;                                      /* service threads that got no core to themselves */
		for (int j = 0; j < n_groups; j++)
			group_free[j] = group_size[j];
		if (bind) {
			int k = emptiest_core(group_cores[0], core_cpus, next_cpu);
			if (k >= 0) {
				thread_cpu[0] = core_cpus[k][next_cpu[k]++];
				group_free[0] -= 1;
			}
		}
		for (int j = 0; j < n_groups; j++) {
			thread_group[1 + j] = j;
			if (not bind)
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
			warnx("MPI: %d dict thread(s) share a core with another service thread: too few cores per"
			      " group for one each", shared);

		/* the emptiest group first, so that no collision queue is left without a producer to drain it */
		for (int s = 0; s < n_producers; s++) {
			int tid = 1 + n_groups + s;
			int best = -1;
			for (int j = 0; bind && j < n_groups; j++)
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

		/* a group whose CPUs ran out borrows a producer from the fullest one: every queue keeps a consumer */
		for (int j = 0; producer_per_dict && j < n_groups; j++) {
			if (n_walk[j] > 0)
				continue;
			int from = 0;
			for (int k = 0; k < n_groups; k++)
				if (n_walk[k] > n_walk[from])
					from = k;
			for (int tid = 1 + n_groups; tid < n_threads; tid++)
				if (thread_group[tid] == from) {
					thread_group[tid] = j;
					n_walk[from] -= 1;
					n_walk[j] += 1;
					break;
				}
			if (rank == 0)
				warnx("MPI: group %d has no CPU left for a producer: one of group %d resolves for it",
				      j, from);
		}
	}

	/* the plan, printed by the startup banner before the team exists */
	void report() const
	{
		printf("MPI: rank %d has %d CPUs in its affinity mask over %d NUMA node(s); threads %s\n",
			rank, n_avail_cpu, n_numa_nodes, bind ? "pinned" : "NOT pinned (--no-bind)");
		int gmin = *std::min_element(group_size.begin(), group_size.end());
		int gmax = *std::max_element(group_size.begin(), group_size.end());
		if (cache_level > 0)
			printf("MPI: %d thread group(s) of %d..%d CPUs over %d L%d domain(s), one dict thread each\n",
				n_groups, gmin, gmax, n_caches, cache_level);
		else
			printf("MPI: %d thread group(s) of %d..%d CPUs; no cache is shared by several cores,"
			       " so the mask is one domain\n", n_groups, gmin, gmax);
		if (bind) {
			printf("MPI: dict threads on CPUs");
			for (int j = 0; j < n_groups; j++)
				printf(" %d", thread_cpu[1 + j]);
			printf("\n");
		}
	}

	/*
	 * What the kernel actually did, once the team is pinned: one line per NUMA node the threads landed
	 * on, then the groups of each cache domain while they are few enough to be worth listing.  `numa`
	 * and `role` are what each thread measured for itself, in tid order (PROTOCOL.md §4.1).
	 */
	void report_measured(const std::vector<int> &numa, const std::vector<int> &role) const
	{
		std::vector<int> numa_ids;
		for (size_t tid = 0; tid < numa.size(); tid++)
			if (std::find(numa_ids.begin(), numa_ids.end(), numa[tid]) == numa_ids.end())
				numa_ids.push_back(numa[tid]);
		std::sort(numa_ids.begin(), numa_ids.end());
		for (size_t k = 0; k < numa_ids.size(); k++) {
			int n_comm = 0, n_dict = 0, n_prod = 0;
			for (size_t tid = 0; tid < numa.size(); tid++) {
				if (numa[tid] != numa_ids[k])
					continue;
				if (role[tid] == COMM)
					n_comm++;
				else if (role[tid] == DICT)
					n_dict++;
				else
					n_prod++;
			}
			printf("NUMA node %d: %s%d dict + %d prod\n", numa_ids[k], n_comm ? "comm + " : "",
				n_dict, n_prod);
		}
		if (cache_level > 0 && n_groups >= n_caches && n_groups <= 32)
			for (int c = 0; c < n_caches; c++) {
				printf("L%d domain %d: groups", cache_level, c);
				for (int j = 0; j < n_groups; j++)
					if (group_cache[j] == c)
						printf(" %d(%d cpu)", j, group_size[j]);
				printf("\n");
			}
		else if (cache_level > 0 && n_groups <= 32)
			for (int j = 0; j < n_groups; j++)
				printf("group %d: %d cpu, L%d domains %d-%d\n", j, group_size[j], cache_level,
					group_cache[j],
					(j + 1 < n_groups) ? group_cache[j + 1] - 1 : n_caches - 1);
		fflush(stdout);
	}
};

}
#endif
