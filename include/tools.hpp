#ifndef MITM_TOOLS
#define MITM_TOOLS

#include <chrono>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <optional>
#include <cerrno>
#include <err.h>
#include <sched.h>
#include <pthread.h>
#include <hwloc.h>

using std::vector;
using std::pair;
using std::tuple;
using std::optional;
using std::nullopt;

#include "types.h"

static_assert(HWLOC_API_VERSION >= 0x00020000, "hwloc 2.x is required (NUMA nodes as memory children)");

namespace mitm {

/* spin-wait hint: these threads are pinned and own their core, so we spin rather
   than yield, but politely. */
static inline void cpu_relax()
{
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#elif defined(__aarch64__)
    __asm__ __volatile__("yield");
#endif
}

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


/* the low n bits; all 64 for n >= 64 */
u64 make_mask(int n)
{
    return (n >= 64) ? 0xffffffffffffffffull : (1ull << n) - 1;
}

/* wall-clock seconds */
double wtime()
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));
  return seconds;
}

/* murmur64 hash, tailored for 64-bit ints.  Cf. Daniel Lemire */
u64 murmur64(u64 h)
{
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    return h;
}

/* one 64-bit hash of a pair, for the HyperLogLog */
u64 murmur128(u64 x, u64 y)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;

    y *= 0xc6a4a7935bd1e995LLU;
    y ^= y >> 47;
    y *= 0xc6a4a7935bd1e995LLU;
    y ^= x;
    y *= 0xc6a4a7935bd1e995LLU;
    return y;
  }


/* deterministic RNG based on TRIVIUM */
class PRNG {
private:
    u64 s11;                    /* TRIVIUM's 93-bit register: low word ... */
    u64 s12;                    /* ... and high word */
    u64 s21;                    /* the 84-bit register: low word ... */
    u64 s22;                    /* ... and high word */
    u64 s31;                    /* the 111-bit register: low word ... */
    u64 s32 = 1;                /* ... and high word */

    void setseed()
    {
        s11 = this->seed;
        s12 = 0;
        s21 = this->seq;
        s22 = 0;
        s31 = 0;
        s32 = 0x700000000000;
        for (int i = 0; i < 18; i++)  /* blank rounds */
            rand();
    }

public:
    const u64 seed;             /* the key */
    const u64 seq;              /* the IV: another stream from the same seed */

    static u64 read_urandom()
    {
        union {
            u64 value;
            char cs[sizeof(u64)];
        } u;
        std::ifstream rfin("/dev/urandom");
        rfin.read(u.cs, sizeof(u.cs));
        rfin.close();
        return u.value;
    }

    u64 rand()
    {
        u64 s66 = (s12 << 62) ^ (s11 >> 2);
        u64 s93 = (s12 << 35) ^ (s11 >> 29);
        u64 s162 = (s22 << 59) ^ (s21 >> 5);
        u64 s177 = (s22 << 44) ^ (s21 >> 20);
        u64 s243 = (s32 << 62) ^ (s31 >> 2);
        u64 s288 = (s32 << 17) ^ (s31 >> 47);
        u64 s91 = (s12 << 37) ^ (s11 >> 27);
        u64 s92 = (s12 << 36) ^ (s11 >> 28);
        u64 s171 = (s22 << 50) ^ (s21 >> 14);
        u64 s175 = (s22 << 46) ^ (s21 >> 18);
        u64 s176 = (s22 << 45) ^ (s21 >> 19);
        u64 s264 = (s32 << 41) ^ (s31 >> 23);
        u64 s286 = (s32 << 19) ^ (s31 >> 45);
        u64 s287 = (s32 << 18) ^ (s31 >> 46);
        u64 s69 = (s12 << 59) ^ (s11 >> 5);
        u64 t1 = s66 ^ s93; /* update */
        u64 t2 = s162 ^ s177;
        u64 t3 = s243 ^ s288;
        u64 z = t1 ^ t2 ^ t3;
        t1 ^= (s91 & s92) ^ s171;
        t2 ^= (s175 & s176) ^ s264;
        t3 ^= (s286 & s287) ^ s69;
        s12 = s11;    /* rotate */
        s11 = t3;
        s22 = s21;
        s21 = t1;
        s32 = s31;
        s31 = t2;
        return z;
    }

    PRNG(u64 seed, u64 seq) : seed(seed), seq(seq)  { setseed(); }
    PRNG(u64 seed) : seed(seed), seq(0) { setseed(); }
    PRNG() : seed(read_urandom()), seq(0) { setseed(); }
};

/* n as a short human string ("1.5G"): at most 7 characters plus the NUL, so an 8-byte target */
void human_format(u64 n, char *target)
{
    if (n < 1000) {
        sprintf(target, "%" PRId64, n);
        return;
    }
    if (n < 1000000) {
        sprintf(target, "%.1fK", n / 1e3);
        return;
    }
    if (n < 1000000000) {
        sprintf(target, "%.1fM", n / 1e6);
        return;
    }
    if (n < 1000000000000ll) {
        sprintf(target, "%.1fG", n / 1e9);
        return;
    }
    if (n < 1000000000000000ll) {
        sprintf(target, "%.1fT", n / 1e12);
        return;
    }
}

/* the inverse of human_format: "4G" -> 4000000000, a bare number as is */
u64 human_parse(const std::string &_h)
{
    std::string h(_h);
    int n = h.length();
    if (h[n - 1] == 'T') {
        h.pop_back();
        return 1000000000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'G') {
        h.pop_back();
        return 1000000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'M') {
        h.pop_back();
        return 1000000ll * (u64) std::stoi(h);
    }
    if (h[n - 1] == 'K') {
        h.pop_back();
        return 1000ll * (u64) std::stoi(h);
    }
    return std::stoull(h);
}
}
#endif
