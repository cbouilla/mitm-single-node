#ifndef MITM_DIRECT_DICT
#define MITM_DIRECT_DICT

#include <cassert>
#include <err.h>

#include "tools.hpp"
#include "direct_common.hpp"

namespace mitm::direct {

/******************************* dictionary shard *****************************/

/*
 * One shard of the distributed dictionary (PROTOCOL.md §8): linear probing over 8-byte slots, insert-only
 * while the round fills, probe-only while it probes, emptied after that.  A slot is the preimage in its
 * low n bits, the key's check bits above, and bit 63 set; zero is empty.  Built by its dict thread once
 * it is pinned (§4.1).
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

	/* where `key`'s run starts.  `key` is the shard's part of the routing key, PROTOCOL.md §2.1 */
	u64 home(u64 key) const
	{
		return key % n_slots;
	}

	/* what every slot of `key`'s entries carries beside its preimage */
	u64 tag(u64 key) const
	{
		return OCCUPIED | (((key / n_slots) & cmask) << n);
	}

	/* the slot after h, around the table */
	u64 step(u64 h) const
	{
		return (h + 1 == n_slots) ? 0 : h + 1;
	}

	/* store (key, x) in the first empty slot of key's run.  `steps` tallies the slots visited */
	void insert(u64 key, u64 x, u64 &steps)
	{
		u64 slot = tag(key) | x;
		u64 h = home(key);
		for (u64 s = 0; s < n_slots; s++) {
			if (A[h] == 0) {
				A[h] = slot;
				steps += s;
				return;
			}
			h = step(h);
		}
		errx(1, "direct: a dictionary shard is full (%" PRIu64 " slots): lower --fill", n_slots);
	}
};


/******************************* the dict thread ******************************/

/*
 * A dict thread's round (PROTOCOL.md §8): while the round fills, every point delivered to its queue is
 * an entry, (key = f(x), x), inserted; while it probes, every point is a probe, (key = g(y), y): each
 * slot of the key's run whose check bits match is a match, verified on the spot with one evaluation of
 * f -- the false positives of the check bits stop there -- and a true match is tested for the golden
 * pair.  Nothing is handed to anybody: a match costs about as much as a probe.  Wind-down: ctx.state,
 * PROTOCOL.md §3.3.
 */
template <class Wrapper>
void Scheme::dict_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
                         SharedContext<Scheme> &shared, int index)
{
	SPSCQueue &in = *ctx.q;
	DirectDict &dict = *shared.shards[index];
	u64 *ctr = ctx.ctr;
	const bool probing = (shared.header.phase == PROBE);
	static const size_t BATCH = 64;
	Point staging[BATCH];

	for (;;) {
		while (not in.empty()) {
			size_t k = in.pop_bulk(staging, BATCH);
			for (size_t t = 0; t < k; t++) {
				const Point &p = staging[t];
				u64 key = p.key / params.n_dicts;
				if (not probing) {
					dict.insert(key, p.val, ctr[N_STEPS]);
					ctr[N_INSERT] += 1;
					continue;
				}
				ctr[N_PROBE] += 1;
				u64 tg = dict.tag(key);
				u64 h = dict.home(key);
				u64 s = 0;
				for (; s < dict.n_slots && dict.A[h] != 0; s++, h = dict.step(h)) {
					if ((dict.A[h] & ~dict.xmask) != tg)
						continue;
					ctr[N_MATCH] += 1;
					u64 x = dict.A[h] & dict.xmask;
					u64 y = p.val;
					if (wrapper.fill(x) != p.key) {
						ctr[BAD_MATCH] += 1;
						continue;
					}
					ctr[N_COLLISIONS] += 1;
					if (not wrapper.good(x, y))
						continue;
					printf("\nFound golden pair! round=%" PRIu64 " x=%" PRIx64 " y=%" PRIx64 "\n",
					       shared.header.round, x, y);
					shared.set_golden(shared.header.round, x, y);
				}
				ctr[N_STEPS] += s;
			}
		}

		/* state first, then emptiness: a point pushed between the two is still consumed (PROTOCOL.md §3.3) */
		if (ctx.state.load(std::memory_order_acquire) == DRAIN) {
			if (in.empty()) {
				ctx.state.store(QUIESCENT, std::memory_order_release);
				return;
			}
			continue;
		}

		cpu_relax();
	}
}

}
#endif
