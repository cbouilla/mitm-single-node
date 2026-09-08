#ifndef MITM_ROUTER_PLACEMENT
#define MITM_ROUTER_PLACEMENT

#include <sched.h>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <err.h>

#include "topology.hpp"

/*
 * Where the Router puts a node's threads, and the groups it forms.  A group is a sub-team of senders and
 * receivers in one cache domain: the staging and the free list will shard on it, and a consumer (PCS) pairs
 * each producer with the dict thread of its group.  By default the Router pins and forms the groups (AUTO);
 * with pin == false it forms them without pinning; with a per-thread color it takes the caller's groups.
 * The role enum and ROUTER_GROUP_AUTO live here, the lowest Router file, so both this and common.hpp see them.
 */

namespace mitm {

enum router_role { ROUTER_SERVICE = 0, ROUTER_SENDER = 1, ROUTER_RECEIVER = 2 };

static constexpr int ROUTER_GROUP_AUTO = -1;   /* Router_Init's `group`: let the Router form the groups */


/*
 * The plan, derived once from the union of the threads' affinity masks and const from then on.  Every
 * per-thread vector is indexed by OpenMP thread id; every per-group vector by group id, domain-major in AUTO.
 */
struct RouterPlacement {
	bool pin = false;                      /* pin the threads and form the groups, or leave the CPUs to the caller */
	int n_threads = 0;                     /* the team: 1 service + S senders + R receivers */
	int n_avail_cpu = 0;                   /* CPUs in the union mask */
	int n_numa_nodes = 0;                  /* NUMA nodes holding one of them */
	int cache_level = 0;                   /* the level a group sits in; 0 == none shared by several cores */
	int n_domains = 0;                     /* cache domains holding a CPU of the mask; 0 when not pinned */
	int n_groups = 0;                      /* G: the groups formed */

	std::vector<int> thread_role;          /* per thread, its router_role: what the reports read */
	std::vector<int> thread_cpu;           /* the CPU to pin it to; -1 == do not pin */
	std::vector<int> thread_domain;        /* the cache domain it lands in; -1 when not pinned */
	std::vector<int> thread_numa;          /* the NUMA node it lands on; -1 when not pinned */
	std::vector<int> thread_group;         /* its group; -1 for the service thread */

	std::vector<int> group_domain;         /* the cache domain of each group; -1 in colors mode */
	std::vector<int> group_nsend;          /* senders in each group */
	std::vector<std::vector<int>> group_receivers;   /* per group, the local indices of its receivers */

	RouterPlacement() = default;

	/*
	 * From the union `mask`, the team's `roles` and `colors` (a worker's group, or ROUTER_GROUP_AUTO), the
	 * options and the rank (for the reports).  Pinned: topology, groups per domain, then round-robin the
	 * receivers and the senders over the domains into the least-filled group.  Unpinned: groups by index.
	 * A color on any worker switches the node to the caller's groups, pinning still the Router's.
	 */
	RouterPlacement(const cpu_set_t &mask, const std::vector<int> &roles, const std::vector<int> &colors,
	                bool pin_, int cache_level_opt, int group_size, int rank)
		: pin(pin_), n_threads((int) roles.size()), n_avail_cpu(CPU_COUNT(&mask)), thread_role(roles)
	{
		thread_cpu.assign(n_threads, -1);
		thread_domain.assign(n_threads, -1);
		thread_numa.assign(n_threads, -1);
		thread_group.assign(n_threads, -1);

		int S = 0;                         /* senders, receivers, and whether every worker asked for AUTO */
		int R = 0;
		bool auto_groups = true;
		for (int t = 0; t < n_threads; t++) {
			if (roles[t] == ROUTER_SENDER)
				S += 1;
			else if (roles[t] == ROUTER_RECEIVER)
				R += 1;
			if (roles[t] != ROUTER_SERVICE && colors[t] != ROUTER_GROUP_AUTO)
				auto_groups = false;
		}

		if (not auto_groups)
			group_by_colors(roles, colors, rank);

		if (not pin) {
			place_unpinned(roles, auto_groups, S, R, group_size);
			return;
		}
		place_pinned(mask, roles, auto_groups, cache_level_opt, group_size, rank);
	}

	/* the group in `groups` with the fewest of `count` so far, ties to the lowest id.  Reached from the
	 * receiver and the sender passes of the pinned AUTO placement. */
	int least_filled(const std::vector<int> &groups, const std::vector<int> &count) const
	{
		int best = groups[0];
		for (size_t k = 1; k < groups.size(); k++)
			if (count[groups[k]] < count[best])
				best = groups[k];
		return best;
	}

	/*
	 * The caller's groups: a worker's group is its color, G the largest plus one, and every id in [0, G) must
	 * be used (contiguous, MPI_Comm_split style).  No cache domain: group_domain is -1.  Reached from the
	 * constructor when any worker passed a color rather than ROUTER_GROUP_AUTO.
	 */
	void group_by_colors(const std::vector<int> &roles, const std::vector<int> &colors, int rank)
	{
		n_groups = 0;
		for (int t = 0; t < n_threads; t++) {
			if (roles[t] == ROUTER_SERVICE)
				continue;
			if (colors[t] < 0)
				errx(1, "Router: rank %d: a group must be set on every worker or none (a worker passed"
				        " ROUTER_GROUP_AUTO while another passed a group)", rank);
			n_groups = std::max(n_groups, colors[t] + 1);
		}
		group_domain.assign(n_groups, -1);
		group_nsend.assign(n_groups, 0);
		group_receivers.assign(n_groups, std::vector<int>());
		std::vector<char> used(n_groups, 0);
		int r_local = 0;
		for (int t = 0; t < n_threads; t++) {
			if (roles[t] == ROUTER_SERVICE)
				continue;
			thread_group[t] = colors[t];
			used[colors[t]] = 1;
			if (roles[t] == ROUTER_SENDER)
				group_nsend[colors[t]] += 1;
			else
				group_receivers[colors[t]].push_back(r_local++);
		}
		for (int g = 0; g < n_groups; g++)
			if (not used[g])
				errx(1, "Router: rank %d: group %d is empty; colors must be 0..G-1 contiguous", rank, g);
	}

	/*
	 * No hwloc: the CPUs stay unset and the groups fall out of the worker index.  AUTO makes
	 * max(ceil(S / group_size), ceil(R / group_size)) groups and slices both roles across them; colors were
	 * already applied.  Reached from the constructor when pin is false.
	 */
	void place_unpinned(const std::vector<int> &roles, bool auto_groups, int S, int R, int group_size)
	{
		n_domains = 0;
		cache_level = 0;
		if (not auto_groups)
			return;
		int gs = (S + group_size - 1) / group_size;
		int gr = (R + group_size - 1) / group_size;
		n_groups = std::max(std::max(gs, gr), 1);
		group_domain.assign(n_groups, -1);
		group_nsend.assign(n_groups, 0);
		group_receivers.assign(n_groups, std::vector<int>());
		int r_local = 0;
		int s_local = 0;
		for (int t = 0; t < n_threads; t++) {
			if (roles[t] == ROUTER_RECEIVER) {
				int g = r_local % n_groups;
				thread_group[t] = g;
				group_receivers[g].push_back(r_local++);
			} else if (roles[t] == ROUTER_SENDER) {
				int g = s_local++ % n_groups;
				thread_group[t] = g;
				group_nsend[g] += 1;
			}
		}
	}

	/*
	 * The pinned placement: the service on domain 0, then the receivers and the senders round-robined over the
	 * domains on ONE cursor that starts after the service's domain, so that the two passes' remainders do not
	 * stack up on the low domains (they did, and the team then doubled up cores at one end of the machine
	 * while leaving cores idle at the other), each onto its domain's emptiest core (a full domain borrows
	 * the globally emptiest core) and, in AUTO, into its domain's least-filled group.  Reached from the
	 * constructor when pin is true.
	 */
	void place_pinned(const cpu_set_t &mask, const std::vector<int> &roles, bool auto_groups, int cache_level_opt,
	                  int group_size, int rank)
	{
		std::vector<int> numa_of_cpu;
		std::vector<int> cache_of_cpu;
		std::vector<int> core_of_cpu;
		cache_level = topology_of_cpus(mask, cache_level_opt, numa_of_cpu, cache_of_cpu, core_of_cpu);

		std::vector<int> numa_ids;                  /* the mask's distinct NUMA nodes, ascending */
		for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
			if (not CPU_ISSET(cpu, &mask))
				continue;
			if (numa_of_cpu[cpu] < 0)
				errx(1, "Router: rank %d: CPU %d is in no NUMA node hwloc knows (use pin = false)", rank, cpu);
			if (std::find(numa_ids.begin(), numa_ids.end(), numa_of_cpu[cpu]) == numa_ids.end())
				numa_ids.push_back(numa_of_cpu[cpu]);
			if (cache_level > 0 && cache_of_cpu[cpu] < 0) {
				if (rank == 0)
					warnx("Router: CPU %d is in no L%d cache hwloc knows: the whole mask becomes one domain",
					      cpu, cache_level);
				cache_level = 0;
			}
		}
		n_numa_nodes = (int) numa_ids.size();

		std::vector<std::vector<int>> core_cpus;    /* the mask's CPUs of each dense core */
		std::vector<std::vector<int>> cache_cores;  /* the dense core indices of each cache domain */
		n_domains = cores_by_domain(mask, cache_level, cache_of_cpu, core_of_cpu, core_cpus, cache_cores);
		std::vector<int> core_domain(core_cpus.size(), 0);   /* each dense core's domain */
		for (int d = 0; d < n_domains; d++)
			for (size_t k = 0; k < cache_cores[d].size(); k++)
				core_domain[cache_cores[d][k]] = d;
		std::vector<int> all_cores(core_cpus.size());        /* every core, for the full-domain fallback */
		for (size_t k = 0; k < core_cpus.size(); k++)
			all_cores[k] = (int) k;

		std::vector<std::vector<int>> domain_groups(n_domains);   /* the group ids of each domain (AUTO) */
		if (auto_groups) {
			n_groups = 0;
			group_domain.clear();
			for (int d = 0; d < n_domains; d++) {
				int cores_d = (int) cache_cores[d].size();
				int i = 1;                  /* groups of at most group_size cores: the smallest i that fits */
				while (cores_d / i > group_size)
					i += 1;
				for (int k = 0; k < i; k++) {
					domain_groups[d].push_back(n_groups);
					group_domain.push_back(d);
					n_groups += 1;
				}
			}
			group_nsend.assign(n_groups, 0);
			group_receivers.assign(n_groups, std::vector<int>());
		}

		std::vector<int> grp_recv(n_groups, 0);     /* receivers placed in each group so far, for least_filled */
		std::vector<size_t> next_cpu(core_cpus.size(), 0);   /* each core's next free CPU */

		int sk = emptiest_core(cache_cores[0], core_cpus, next_cpu);   /* the service: domain 0, a core of its own */
		thread_cpu[0] = core_cpus[sk][next_cpu[sk]++];
		thread_domain[0] = 0;
		thread_numa[0] = numa_of_cpu[thread_cpu[0]];

		int dc = n_domains > 1 ? 1 : 0;             /* the domain cursor: the service already took a core of domain 0 */
		int r_local = 0;
		for (int t = 0; t < n_threads; t++) {
			if (roles[t] != ROUTER_RECEIVER)
				continue;
			int d = dc;
			dc = (dc + 1) % n_domains;
			int k = emptiest_core(cache_cores[d], core_cpus, next_cpu);
			if (k < 0) {                    /* domain d is full: take the globally emptiest core */
				k = emptiest_core(all_cores, core_cpus, next_cpu);
				d = core_domain[k];
			}
			thread_cpu[t] = core_cpus[k][next_cpu[k]++];
			thread_domain[t] = d;
			thread_numa[t] = numa_of_cpu[thread_cpu[t]];
			if (auto_groups) {
				int g = least_filled(domain_groups[d], grp_recv);
				thread_group[t] = g;
				grp_recv[g] += 1;
				group_receivers[g].push_back(r_local);
			}
			r_local += 1;
		}

		for (int t = 0; t < n_threads; t++) {   /* dc carries on from the receivers: one cursor over the whole
		                                         * team, so the two passes' remainders do not stack up on the low
		                                         * domains, on top of the service's core */
			if (roles[t] != ROUTER_SENDER)
				continue;
			int d = dc;
			dc = (dc + 1) % n_domains;
			int k = emptiest_core(cache_cores[d], core_cpus, next_cpu);
			if (k < 0) {
				k = emptiest_core(all_cores, core_cpus, next_cpu);
				d = core_domain[k];
			}
			thread_cpu[t] = core_cpus[k][next_cpu[k]++];
			thread_domain[t] = d;
			thread_numa[t] = numa_of_cpu[thread_cpu[t]];
			if (auto_groups) {
				int g = least_filled(domain_groups[d], group_nsend);
				thread_group[t] = g;
				group_nsend[g] += 1;
			}
		}
	}

	/* the intended layout, once.  Reached from connect() before the sizes banner, on rank 0 when verbose. */
	void report() const
	{
		if (not pin) {
			printf("Router: threads NOT pinned (pin = false); %d group(s) by worker index\n", n_groups);
			return;
		}
		printf("Router: %d CPUs over %d NUMA node(s)", n_avail_cpu, n_numa_nodes);
		if (cache_level > 0)
			printf(", %d L%d domain(s)", n_domains, cache_level);
		else
			printf(", no cache shared by several cores (one domain)");
		int smin = group_nsend.empty() ? 0 : group_nsend[0];
		int smax = 0;
		int rmin = group_receivers.empty() ? 0 : (int) group_receivers[0].size();
		int rmax = 0;
		for (int g = 0; g < n_groups; g++) {
			smin = std::min(smin, group_nsend[g]);
			smax = std::max(smax, group_nsend[g]);
			rmin = std::min(rmin, (int) group_receivers[g].size());
			rmax = std::max(rmax, (int) group_receivers[g].size());
		}
		printf("; %d group(s) of %d..%d senders and %d..%d receivers, threads pinned\n",
		       n_groups, smin, smax, rmin, rmax);
		if (n_numa_nodes != 1) {
			printf("***** WARNING *****\n");
			printf("---> a rank spans %d NUMA nodes: one MPI rank per NUMA node is wanted\n", n_numa_nodes);
			printf("---> (mpirun --map-by numa --bind-to numa, or the launcher's equivalent)\n");
			printf("***** WARNING *****\n");
		}
	}

	/*
	 * What the kernel honoured, once the team is pinned and every thread's getcpu has confirmed it: per NUMA
	 * node the counts by role, then per domain its groups with their sender/receiver counts while the domains
	 * are few enough to list.  Reached from thread 0's Router_thread constructor, rank 0 when verbose.
	 */
	void report_measured() const
	{
		if (not pin)
			return;
		std::vector<int> ids;
		for (int t = 0; t < n_threads; t++)
			if (std::find(ids.begin(), ids.end(), thread_numa[t]) == ids.end())
				ids.push_back(thread_numa[t]);
		std::sort(ids.begin(), ids.end());
		for (size_t j = 0; j < ids.size(); j++) {
			int n_comm = 0;
			int n_send = 0;
			int n_recv = 0;
			for (int t = 0; t < n_threads; t++) {
				if (thread_numa[t] != ids[j])
					continue;
				if (thread_role[t] == ROUTER_SERVICE)
					n_comm += 1;
				else if (thread_role[t] == ROUTER_SENDER)
					n_send += 1;
				else
					n_recv += 1;
			}
			printf("NUMA node %d: %s%d senders + %d receivers\n", ids[j], n_comm ? "service + " : "",
			       n_send, n_recv);
		}
		if (cache_level > 0 && n_domains <= 32)
			for (int d = 0; d < n_domains; d++) {
				int cnt = 0;                /* colors' groups have no domain (group_domain -1): nothing to list */
				for (int g = 0; g < n_groups; g++)
					cnt += group_domain[g] == d;
				if (cnt == 0)
					continue;
				printf("L%d domain %d: groups", cache_level, d);
				for (int g = 0; g < n_groups; g++)
					if (group_domain[g] == d)
						printf(" %d(%ds %zur)", g, group_nsend[g], group_receivers[g].size());
				printf("\n");
			}
		fflush(stdout);
	}
};

}
#endif
