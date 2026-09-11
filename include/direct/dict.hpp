#ifndef MITM_DIRECT_DICT
#define MITM_DIRECT_DICT

#include <cassert>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cinttypes>
#include <algorithm>
#include <sys/mman.h>

#include "tools.hpp"
#include "router/router.hpp"
#include "direct/params.hpp"
#include "direct/shared.hpp"

/* One shard of the distributed dictionary of preimages, and the thread that owns it. */

namespace mitm::direct {

/******************************* the shard's pages ****************************/

/* the page these mappings ask for: the kernel rounds a MAP_HUGETLB length up to it, munmap does not */
static constexpr size_t HUGE_PAGE = 2 * 1024 * 1024;

/*
 * Map `nbytes` of zeroed memory for a shard: 2 MB pages when the kernel has some reserved, ordinary ones
 * otherwise, `why` saying why MAP_HUGETLB was refused, 0 if it was not.  No MAP_LOCKED: the kernel clears
 * VM_LOCKED on a hugetlb mapping, whose pages cannot be swapped anyway, so the flag locks nothing and its
 * RLIMIT_MEMLOCK check only refuses shards larger than the limit.
 */
static u64 *map_shard(size_t nbytes, int *why)
{
	int flags = MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB | MAP_HUGE_2MB;
	void *ptr = mmap(NULL, nbytes, PROT_READ | PROT_WRITE, flags, -1, 0);
	if (ptr != MAP_FAILED)
		return (u64 *) ptr;
	*why = errno;
	ptr = mmap(NULL, nbytes, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
	if (ptr == MAP_FAILED)
		err(1, "direct: cannot map %zu bytes for a dictionary shard", nbytes);
	madvise(ptr, nbytes, MADV_HUGEPAGE);   /* transparent huge pages, all the kernel has left to offer */
	return (u64 *) ptr;
}

/*
 * How much of [ptr, ptr + nbytes) the kernel really backs with 2 MB pages, read from /proc/self/smaps once
 * every page has been touched.  Exact for MAP_HUGETLB, whose mapping stands alone; an upper bound for
 * transparent pages, which the kernel reports over a region it may have merged with the neighbouring ones.
 */
static u64 huge_bytes(const void *ptr, size_t nbytes)
{
	FILE *f = fopen("/proc/self/smaps", "r");
	if (f == NULL)
		return 0;
	uintptr_t base = (uintptr_t) ptr;
	uintptr_t end = base + nbytes;
	u64 overlap = 0;          /* how much of the region being read falls inside the shard */
	u64 total = 0;            /* the bytes on 2 MB pages so far */
	char line[256];
	while (fgets(line, sizeof(line), f) != NULL) {
		uintptr_t lo, hi;
		int kb;
		if (sscanf(line, "%" SCNxPTR "-%" SCNxPTR, &lo, &hi) == 2) {
			overlap = (lo < end && base < hi) ? std::min(end, hi) - std::max(base, lo) : 0;
			continue;
		}
		if (overlap == 0)
			continue;
		if (sscanf(line, "AnonHugePages: %d kB", &kb) == 1)
			total += (u64) kb * 1024;
		if (sscanf(line, "MMUPageSize: %d kB", &kb) == 1 && kb >= 2048)
			total += overlap;      /* a hugetlb region: every page of it is a huge one */
	}
	fclose(f);
	return std::min(total, (u64) nbytes);
}


/******************************* dictionary shard *****************************/

/*
 * One shard of the distributed dictionary: linear probing over 8-byte slots, insert-only while the round fills,
 * probe-only while it probes, emptied after that.  A slot is the preimage in its low n bits, the key's check
 * bits above, and bit 63 set; zero is empty.  The key is the hashed image the Router delivered, mixed once more
 * here: the routing consumed the top bits of its low word, and a shard above 2^32 slots cut from the key itself
 * would leave part of its slots unreachable.  Built by its dict thread once it is pinned.
 */
class DirectDict {
public:
	static constexpr u64 OCCUPIED = 1ull << 63;   /* the occupancy bit, above the check bits */

	const u64 n_slots;     /* size of A */
	const int n;           /* preimage bits */
	const u64 xmask;       /* the low n bits: the preimage */
	const u64 cmask;       /* the check bits, 63 - n of them, before their shift */
	u64 *A;                /* A[i] == OCCUPIED | check << n | preimage, or 0 */
	size_t nbytes;         /* the mapping: the slots rounded up to a whole 2 MB page, so up to that much over */
	u64 huge;              /* of those bytes, what the kernel really backs with 2 MB pages */
	int hugetlb_errno;     /* why MAP_HUGETLB was refused; 0 if it was not */

	DirectDict(u64 n_slots, int n) : n_slots(n_slots), n(n), xmask(make_mask(n)), cmask(make_mask(63 - n)),
	                                 A(nullptr), nbytes(0), huge(0), hugetlb_errno(0)
	{
		assert(n <= 63);
		if (n_slots == 0)
			return;            /* the service thread and the producers hold an empty shard */
		nbytes = (n_slots * sizeof(u64) + HUGE_PAGE - 1) & ~(HUGE_PAGE - 1);
		A = map_shard(nbytes, &hugetlb_errno);
		memset(A, 0, nbytes);      /* by the calling thread: the first touch of every page, and its promotion */
		huge = huge_bytes(A, nbytes);
	}

	~DirectDict()
	{
		if (A != nullptr)
			munmap(A, nbytes);
	}

	DirectDict(const DirectDict &) = delete;
	DirectDict &operator=(const DirectDict &) = delete;

	/*
	 * Store (key, x) in the first empty slot of its run, which starts at the top bits of the mixed key, scaled
	 * to the table by one multiplication.  Tallies the insert and the slots it visited.
	 */
	void insert(u64 key, u64 x, u64 *ctr)
	{
		u64 h = murmur64(key);
		u64 slot = OCCUPIED | ((h & cmask) << n) | x;
		u64 s = (u64) (((unsigned __int128) h * n_slots) >> 64);
		for (u64 k = 0; k < n_slots; k++, s = (s + 1 == n_slots) ? 0 : s + 1) {
			if (A[s] != 0)
				continue;
			A[s] = slot;
			ctr[N_INSERT] += 1;
			ctr[N_STEPS] += k;
			return;
		}
		errx(1, "direct: a dictionary shard is full (%" PRIu64 " slots): lower --fill", n_slots);
	}

	/*
	 * One probe (key, y) against the shard: every slot of the key's run carrying its tag -- bit 63 and the
	 * mixed key's low check bits -- is a match, kept only if the preimage really maps to the key, and a true
	 * collision is tested for the golden pair, which goes to `shared`.
	 */
	template <class Wrapper>
	void probe(const Wrapper &wrapper, Shared &shared, u64 *ctr, u64 round, u64 key, u64 y) const
	{
		u64 h = murmur64(key);
		u64 tag = OCCUPIED | ((h & cmask) << n);
		u64 s = (u64) (((unsigned __int128) h * n_slots) >> 64);
		u64 k = 0;
		for (; k < n_slots && A[s] != 0; k++, s = (s + 1 == n_slots) ? 0 : s + 1) {
			if ((A[s] & ~xmask) != tag)
				continue;
			ctr[N_MATCH] += 1;
			u64 x = A[s] & xmask;
			if (murmur64(wrapper.pb.f(x)) != key) {
				ctr[BAD_MATCH] += 1;
				continue;
			}
			ctr[N_COLLISIONS] += 1;
			if (not wrapper.good(x, y))
				continue;
			printf("\nFound golden pair! round=%" PRIu64 " x=%" PRIx64 " y=%" PRIx64 "\n", round, x, y);
			shared.set_golden(x, y);
		}
		ctr[N_PROBE] += 1;
		ctr[N_STEPS] += k;
	}

	/* empty the shard after a PROBE phase; the owner alone today */
	void flush()
	{
		memset(A, 0, n_slots * sizeof(u64));
	}
};


/******************************* the dict thread ******************************/

/*
 * A dict thread's phase: take the points the Router delivers, one at a time, until nothing more can arrive.
 * While the round fills, a point is an entry, (hashed image, preimage), inserted; while it probes, a point is a
 * probe, resolved on the spot: a match costs about what a probe does, so nothing is handed to anybody.  After a
 * PROBE phase the shard is emptied.
 */
template <class Wrapper>
void dict_round(Router_thread &rt, const Wrapper &wrapper, Shared &shared, DirectDict &dict, u64 *ctr,
                u64 round, int phase)
{
	for (;;) {
		u64 key;                       /* the hashed image the producer routed on */
		u64 val;                       /* its preimage */
		if (not Router_Pop(&key, &val, rt)) {
			if (Router_Test_drained(rt))
				break;
			cpu_relax();               /* no CAS in Router_Pop --> no need for backoff */
			continue;
		}
		if (phase == FILL)
			dict.insert(key, val, ctr);
		else
			dict.probe(wrapper, shared, ctr, round, key, val);
	}
	if (phase == PROBE)
		dict.flush();
}

}
#endif
