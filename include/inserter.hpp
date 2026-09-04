#ifndef MITM_INSERTER
#define MITM_INSERTER

#include <cassert>

#include "tools.hpp"
#include "parameters.hpp"
#include "comm.hpp"

namespace mitm {

/******************************* dictionary shard *****************************/

/*
 * One shard of the distributed dictionary: probed with a trail's endpoint, it returns the chain index
 * (and the length, if it fit) of an earlier trail with the same endpoint.  Built by its inserter once
 * it is pinned (PROTOCOL.md §4.1).
 */
class PcsDict {
public:
	u64 jbits;             /* low bits of a slot: the chain index */
	u64 lbits;             /* jbits + 8: the trail length sits above the chain index */
	u64 jmask;             /* the low jbits */
	u64 lmask;             /* 0xff: lengths saturate at 255 */
	u64 key_mask;          /* the bits above lbits */
	const u64 n_slots;     /* size of A */

	vector<u64> A;         /* A[i][0:jbits] == j.  A[i][jbits:lbits] == len1.  A[lbits:64] == key bits */

	PcsDict(u64 jbits, u64 w) : jbits(jbits), n_slots(w)
	{
		assert(jbits <= 56);
		jmask = make_mask(jbits);
		lmask = make_mask(8);
		lbits = jbits + 8;
		key_mask = (lbits == 64) ? 0 : 0xffffffffffffffff << lbits;
		A.resize(n_slots);     /* zero-filled by the calling thread: the first touch of every page */
	}

	/* zero the shard between rounds; the owner alone today (a collective flush is future work: CLAUDE.md) */
	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i] = 0;
	}

	/*
	 * Probe-and-insert: the (chain index, length) already at this endpoint's slot, if its key matches;
	 * length 0 == it had saturated.  A slot is overwritten only by a trail at least as long.
	 */
	optional<pair<u64, u64>> pop_insert(u64 end, u64 start, u64 len0)
	{
		u64 idx = end % n_slots;
		u64 key = (end / n_slots) << lbits;

		u64 e = A[idx];
		u64 ekey = e & key_mask;
		u64 elen = (e >> jbits) & lmask;

		if (e == 0 || len0 >= elen) {
			if (len0 > lmask)
				len0 = lmask;
			A[idx] = start ^ (len0 << jbits) ^ key;
		}

		if (ekey != key || e == 0)
			return nullopt;

		if (elen == lmask)
			elen = 0;

		return optional(pair(e & jmask, elen));
	}
};

/******************************* the inserter thread **************************/

/*
 * An inserter thread: probes every DP delivered to its queue into its shard and hands each hit to the
 * walkers as a CollisionCandidate.  Wind-down: ctx.state, PROTOCOL.md §3.3.
 *
 * Not done: software-pipelining the probes (prefetch a fresh point, probe one queued a dozen earlier,
 * so that many misses are outstanding at once).  Measured at about 5% on a laptop with a 256 MB
 * dictionary, not worth its ring buffer; retry on cluster hardware before dismissing it.
 */
inline void inserter_thread(ThreadContext &ctx, const Parameters &params, SharedContext &shared,
                            int inserter_index)
{
	SPSCQueue &in = *ctx.q;
	PcsDict &dict = *shared.shards[inserter_index];
	CollisionQueue &coll_q = shared.coll_q;
	u64 *ctr = ctx.ctr;
	static const size_t BATCH = 64;
	DP staging[BATCH];

	for (;;) {
		while (not in.empty()) {
			size_t k = in.pop_bulk(staging, BATCH);
			for (size_t t = 0; t < k; t++) {
				const DP &p = staging[t];
				u64 key = p.x / params.n_inserters;
				ctr[N_PROBE] += 1;
				auto hit = dict.pop_insert(key, p.seed, p.len);
				if (not hit) {
					ctr[BAD_PROBE] += 1;
					continue;
				}
				auto [seed1, len1_maybe] = *hit;
				CollisionCandidate c = {shared.i, p.seed, p.len, key, seed1, len1_maybe};
				if (not coll_q.push(c))
					ctr[DROP_COLL] += 1;
			}
		}

		/* state first, then emptiness: a point pushed between the two is still probed (PROTOCOL.md §3.3) */
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
