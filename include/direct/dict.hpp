#ifndef MITM_DIRECT_DICT
#define MITM_DIRECT_DICT

#include <cassert>
#include <cstring>
#include <cstdio>

#include "tools.hpp"
#include "router/router.hpp"
#include "direct/params.hpp"
#include "direct/shared.hpp"

/* One shard of the distributed dictionary of preimages, and the thread that owns it. */

namespace mitm::direct {

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
	vector<u64> A;         /* A[i] == OCCUPIED | check << n | preimage, or 0 */

	DirectDict(u64 n_slots, int n) : n_slots(n_slots), n(n), xmask(make_mask(n)), cmask(make_mask(63 - n))
	{
		assert(n <= 63);
		A.resize(n_slots);     /* zero-filled by the calling thread: the first touch of every page */
	}

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
		memset(A.data(), 0, n_slots * sizeof(u64));
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
