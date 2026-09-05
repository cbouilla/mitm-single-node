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
	 * length 0 == it had saturated.  A slot is overwritten only by a trail at least as long.  `len0`
	 * arrives saturated by the wire format already, so the clamp below only narrows 8 bits further.
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
 * An inserter thread: probes every DP delivered to its queue into its shard and hands the hits to the
 * walkers of its own group, over its own collision queue and in runs (PROTOCOL.md §3.2).  A run is
 * handed over as soon as there is nothing left to probe, so a partial one never waits on the next hit.
 * Wind-down: ctx.state, PROTOCOL.md §3.3.
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
	CollisionQueue &coll_q = *shared.coll_q[inserter_index];
	u64 *ctr = ctx.ctr;
	u64 jmask = make_mask(params.jbits);
	static const size_t BATCH = 64;
	DP staging[BATCH];
	CollisionCandidate pending[BATCH];    /* hits waiting for the run to be handed over */
	size_t n_pending = 0;

	for (;;) {
		while (not in.empty()) {
			size_t k = in.pop_bulk(staging, BATCH);
			for (size_t t = 0; t < k; t++) {
				const DP &p = staging[t];
				u64 key = p.x / params.n_inserters;
				u64 seed0 = p.jl & jmask;
				u64 len = (p.jl >> params.jbits) & params.len_sat;
				u64 len0_maybe = (len == params.len_sat) ? 0 : len;
				/* unknown on the wire stays unknown in the slot, whatever the 8-bit field could hold */
				u64 dict_len = len0_maybe ? len0_maybe : 0xffffffffffffffffull;
				ctr[N_PROBE] += 1;
				auto hit = dict.pop_insert(key, seed0, dict_len);
				if (not hit) {
					ctr[BAD_PROBE] += 1;
					continue;
				}
				auto [seed1, len1_maybe] = *hit;
				pending[n_pending] = {shared.i, seed0, len0_maybe, p.x, seed1, len1_maybe};
				n_pending += 1;
				if (n_pending == BATCH) {
					ctr[DROP_COLL] += n_pending - coll_q.push_bulk(pending, n_pending);
					n_pending = 0;
				}
			}
		}

		/* nothing left to probe: the run goes over now, so that the drain never has to flush one (§4.5) */
		if (n_pending > 0) {
			ctr[DROP_COLL] += n_pending - coll_q.push_bulk(pending, n_pending);
			n_pending = 0;
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
