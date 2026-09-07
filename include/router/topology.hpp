#ifndef MITM_ROUTER_TOPOLOGY
#define MITM_ROUTER_TOPOLOGY

#include <err.h>
#include <sched.h>
#include <pthread.h>
#include <cerrno>
#include <algorithm>
#include <vector>
#include <hwloc.h>

/*
 * The hwloc primitives shared by the engine's placement and the Router's: pinning, the cache level a
 * thread group sits in, the NUMA node / cache domain / core of every CPU of a mask, the mask's cores by
 * domain, and the emptiest core of a set.  What a CPU is, and nothing about who runs on it.
 */

static_assert(HWLOC_API_VERSION >= 0x00020000, "hwloc 2.x is required (NUMA nodes as memory children)");

namespace mitm {

/* pin the calling thread to `cpu`.  Returns cpu, or -1 with errno set: a failure, or cpu < 0 */
static inline int pin_to_cpu(int cpu)
{
	if (cpu < 0)
		return -1;
	cpu_set_t set;
	CPU_ZERO(&set);
	CPU_SET(cpu, &set);
	int rc = pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
	if (rc != 0) {
		errno = rc;
		return -1;
	}
	return cpu;
}

/*
 * The lowest cache level shared by several CORES, which is what a thread group is built around; 0 when
 * every cache is private to one core.  Depths run from the machine down to the PUs, so walking them
 * backwards visits the caches from the smallest up and stops at the first shared one.
 */
static inline int lowest_shared_cache_level(hwloc_topology_t topo)
{
	for (int depth = hwloc_topology_get_depth(topo) - 1; depth >= 0; depth--) {
		hwloc_obj_t obj = hwloc_get_obj_by_depth(topo, depth, 0);
		if (obj == NULL || not hwloc_obj_type_is_dcache(obj->type))
			continue;
		if (hwloc_get_nbobjs_inside_cpuset_by_type(topo, obj->cpuset, HWLOC_OBJ_CORE) > 1)
			return (int) obj->attr->cache.depth;
	}
	return 0;
}

/*
 * NUMA node (hwloc os_index), cache domain and physical core (both a dense index, ascending) of every
 * CPU of `mask`, from one topology load; -1 outside the mask, or in no such object.  `level` is the
 * cache level that defines a domain, or 0 to take the lowest one shared by several cores; the level
 * used is returned.  The core map is what makes SMT siblings visible to the placement.
 * No topology flag on purpose: RESTRICT_TO_CPUBINDING would defeat the HWLOC_SYNTHETIC test path.
 */
static inline int topology_of_cpus(const cpu_set_t &mask, int level, std::vector<int> &numa_node_of_cpu,
                                   std::vector<int> &cache_of_cpu, std::vector<int> &core_of_cpu)
{
	numa_node_of_cpu.assign(CPU_SETSIZE, -1);
	cache_of_cpu.assign(CPU_SETSIZE, -1);
	core_of_cpu.assign(CPU_SETSIZE, -1);
	hwloc_topology_t topo;
	if (hwloc_topology_init(&topo) != 0)
		err(1, "hwloc_topology_init");
	if (hwloc_topology_load(topo) != 0)
		err(1, "hwloc_topology_load");

	hwloc_obj_t node = NULL;
	while ((node = hwloc_get_next_obj_by_type(topo, HWLOC_OBJ_NUMANODE, node)) != NULL) {
		unsigned cpu;
		hwloc_bitmap_foreach_begin(cpu, node->cpuset)
			if (cpu < CPU_SETSIZE && CPU_ISSET(cpu, &mask) && numa_node_of_cpu[cpu] < 0)
				numa_node_of_cpu[cpu] = (int) node->os_index;
		hwloc_bitmap_foreach_end();
	}

	hwloc_obj_t core = NULL;
	int n_core = 0;
	while ((core = hwloc_get_next_obj_by_type(topo, HWLOC_OBJ_CORE, core)) != NULL) {
		bool used = false;
		unsigned cpu;
		hwloc_bitmap_foreach_begin(cpu, core->cpuset)
			if (cpu < CPU_SETSIZE && CPU_ISSET(cpu, &mask) && core_of_cpu[cpu] < 0) {
				core_of_cpu[cpu] = n_core;
				used = true;
			}
		hwloc_bitmap_foreach_end();
		if (used)
			n_core += 1;
	}

	if (level == 0)
		level = lowest_shared_cache_level(topo);
	int id = 0;                        /* only domains holding a CPU of the mask are numbered */
	for (int depth = hwloc_topology_get_depth(topo) - 1; level > 0 && depth >= 0; depth--) {
		unsigned nobj = hwloc_get_nbobjs_by_depth(topo, depth);
		for (unsigned k = 0; k < nobj; k++) {
			hwloc_obj_t obj = hwloc_get_obj_by_depth(topo, depth, k);
			if (not hwloc_obj_type_is_dcache(obj->type) || (int) obj->attr->cache.depth != level)
				continue;
			bool used = false;
			unsigned cpu;
			hwloc_bitmap_foreach_begin(cpu, obj->cpuset)
				if (cpu < CPU_SETSIZE && CPU_ISSET(cpu, &mask) && cache_of_cpu[cpu] < 0) {
					cache_of_cpu[cpu] = id;
					used = true;
				}
			hwloc_bitmap_foreach_end();
			if (used)
				id += 1;
		}
	}
	hwloc_topology_destroy(topo);
	return (id > 0) ? level : 0;
}

/*
 * The mask's CPUs by dense core index, and the dense core indices by cache domain; one domain (0) when no
 * cache is shared.  Returns the number of domains.  `cache_of_cpu` and `core_of_cpu` are as topology_of_cpus
 * filled them; the per-core domain is used only to bin the cores and is not returned.
 */
static inline int cores_by_domain(const cpu_set_t &mask, int cache_level, const std::vector<int> &cache_of_cpu,
                                  const std::vector<int> &core_of_cpu,
                                  std::vector<std::vector<int>> &core_cpus,
                                  std::vector<std::vector<int>> &cache_cores)
{
	core_cpus.clear();
	std::vector<int> core_cache;                /* cache domain of each core */
	std::vector<int> dense(CPU_SETSIZE, -1);    /* hwloc's core index -> ours */
	int n_caches = 0;
	for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
		if (!CPU_ISSET(cpu, &mask))
			continue;
		int c = (cache_level > 0) ? cache_of_cpu[cpu] : 0;
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
	cache_cores.assign(n_caches, std::vector<int>());
	for (size_t k = 0; k < core_cpus.size(); k++)
		cache_cores[core_cache[k]].push_back(k);
	return n_caches;
}

/*
 * The group's emptiest core that still has a free CPU, or -1.  Handing every thread the emptiest core
 * is what puts the comm thread and the dict threads on cores of their own, and then fills the siblings
 * they left with producers: a producer is compute-bound and fills the issue slots a thread waiting on
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

}
#endif
