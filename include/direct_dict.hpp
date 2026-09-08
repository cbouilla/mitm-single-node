#ifndef MITM_DIRECT_DICT
#define MITM_DIRECT_DICT

#include <cassert>
#include <err.h>

#include "tools.hpp"
#include "direct_common.hpp"

namespace mitm::direct {

/******************************* dictionary shard *****************************/

/*
 * One shard of the distributed dictionary: linear probing over 8-byte slots, insert-only while the round
 * fills, probe-only while it probes, emptied after that.  A slot is the preimage in its low n bits, the key's
 * check bits above, and bit 63 set; zero is empty.  The key is the hashed image the Router delivered, mixed
 * once more here: the routing consumed the top bits of its low word, and a shard above 2^32 slots cut from the
 * key itself would leave part of its slots unreachable.  Built by its dict thread once it is pinned.
 */
class DirectDict {
public:
	static constexpr u64 OCCUPIED = 1ull << 63;

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

	/* zero the shard after a PROBE phase; the owner alone today */
	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i] = 0;
	}

	/* where the run of mixed key h starts: its top bits, scaled to the table by one multiplication */
	u64 home(u64 h) const
	{
		return (u64) (((unsigned __int128) h * n_slots) >> 64);
	}

	/* what every slot holding mixed key h carries beside its preimage: bit 63 and h's low check bits */
	u64 tag(u64 h) const
	{
		return OCCUPIED | ((h & cmask) << n);
	}

	/* the slot after s, around the table */
	u64 step(u64 s) const
	{
		return (s + 1 == n_slots) ? 0 : s + 1;
	}

	/* store (key, x) in the first empty slot of key's run.  `steps` tallies the slots visited */
	void insert(u64 key, u64 x, u64 &steps)
	{
		u64 h = murmur64(key);
		u64 slot = tag(h) | x;
		u64 s = home(h);
		for (u64 k = 0; k < n_slots; k++) {
			if (A[s] == 0) {
				A[s] = slot;
				steps += k;
				return;
			}
			s = step(s);
		}
		errx(1, "direct: a dictionary shard is full (%" PRIu64 " slots): lower --fill", n_slots);
	}
};


/******************************* the dict thread ******************************/

/*
 * One probe (key, y) against the shard: every slot of key's run carrying its tag is a match, kept only if the
 * preimage really maps to the key, and a true collision is tested for the golden pair, which goes to `shared`.
 */
template <class Wrapper>
static void probe(const Wrapper &wrapper, Shared &shared, DirectDict &dict, u64 *ctr, u64 round, u64 key, u64 y)
{
	u64 h = murmur64(key);
	u64 tg = dict.tag(h);
	u64 s = dict.home(h);
	u64 k = 0;
	for (; k < dict.n_slots && dict.A[s] != 0; k++, s = dict.step(s)) {
		if ((dict.A[s] & ~dict.xmask) != tg)
			continue;
		ctr[N_MATCH] += 1;
		u64 x = dict.A[s] & dict.xmask;
		if (murmur64(wrapper.fill(x)) != key) {
			ctr[BAD_MATCH] += 1;
			continue;
		}
		ctr[N_COLLISIONS] += 1;
		u64 yy = y;                    /* good() may swap the pair into the order the problem wants */
		if (not wrapper.good(x, yy))
			continue;
		printf("\nFound golden pair! round=%" PRIu64 " x=%" PRIx64 " y=%" PRIx64 "\n", round, x, yy);
		shared.set_golden(round, x, yy);
	}
	ctr[N_STEPS] += k;
}

/*
 * A dict thread's phase: every block the Router delivers is read where it lies.  While the round fills, a
 * point is an entry, (hashed image, preimage), inserted; while it probes, a point is a probe, resolved on the
 * spot: a match costs about what a probe does, so nothing is handed to anybody.  After a PROBE phase the
 * shard is emptied.
 */
template <class Wrapper>
void dict_round(Router_thread &rt, const Wrapper &wrapper, Shared &shared, DirectDict &dict, u64 *ctr, u64 round,
                int phase)
{
	for (;;) {
		const u64 *pts;
		size_t k = Router_Grab(&pts, rt);
		if (k == 0) {
			if (Router_Test_drained(rt))
				break;
			cpu_relax();
			continue;
		}
		if (phase == FILL) {
			for (size_t i = 0; i < k; i++)
				dict.insert(pts[2 * i], pts[2 * i + 1], ctr[N_STEPS]);
			ctr[N_INSERT] += k;
		} else {
			for (size_t i = 0; i < k; i++)
				probe(wrapper, shared, dict, ctr, round, pts[2 * i], pts[2 * i + 1]);
			ctr[N_PROBE] += k;
		}
		Router_Release(rt);
	}
	if (phase == PROBE)
		dict.flush();
}

}
#endif
