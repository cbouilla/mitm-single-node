#ifndef MITM_DIRECT_PRODUCER
#define MITM_DIRECT_PRODUCER

#include <cassert>

#include "tools.hpp"
#include "direct_common.hpp"

namespace mitm::direct {

/*
 * A producer's round (PROTOCOL.md §8): evaluate the phase's function on its own contiguous sub-range of
 * the phase's domain -- the round's chunk while filling, the whole domain while probing -- vlen inputs
 * at a time, and ship every image as a point, (image, preimage), waiting for room rather than dropping
 * (§5).  Once the range is exhausted it holds itself (RUNNING -> HELD, §3.3) and waits to be drained.
 * The sub-ranges of the n_producers producers cut the phase's domain into equal contiguous pieces, so
 * the partition is the same on every node and needs no message.
 */
template <class Wrapper>
void Scheme::producer_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
                             SharedContext<Scheme> &shared, int index)
{
	constexpr int vlen = Wrapper::vlen;
	SPSCQueue &out = *ctx.q;
	u64 *ctr = ctx.ctr;
	const u64 phase = shared.header.phase;

	/* the phase's domain: chunk `round` of the preimages while filling, all of them while probing */
	u64 lo = 0;
	u64 hi = params.domain;
	if (phase == FILL) {
		lo = shared.header.round * params.per_round;
		hi = std::min(lo + params.per_round, params.domain);
	}

	/* this producer's piece: the global index cuts [lo, hi) into n_producers pieces, the first ones one longer */
	u64 p = (u64) params.rank * params.producers_per_node + index;
	u64 span = hi - lo;
	u64 piece = span / params.n_producers;
	u64 extra = span % params.n_producers;
	u64 my_lo = lo + piece * p + std::min(p, extra);
	u64 my_hi = my_lo + piece + ((p < extra) ? 1 : 0);

	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	size_t since_check = 0;             /* vectors since the state was last read */

	for (u64 base = my_lo; base < my_hi; base += vlen) {
		if (since_check == params.chunk_size) {
			since_check = 0;
			int st = ctx.state.load(std::memory_order_acquire);
			if (st != RUNNING)
				break;                  /* HOLD: the round is over before the range is (§4.5, step 1) */
		}
		since_check += 1;

		int valid = (my_hi - base < (u64) vlen) ? (int) (my_hi - base) : vlen;
		for (int k = 0; k < vlen; k++)
			x[k] = (k < valid) ? base + k : base;   /* a lane past the end still gets an input in the domain */
		wrapper.veval(phase, x, y);
		ctr[N_EVAL] += valid;

		for (int k = 0; k < valid; k++) {
			Point pt = {y[k], x[k]};
			while (not out.push(pt))
				cpu_relax();            /* lossless: wait for the comm thread, which never waits for us (§5) */
		}
		ctr[N_POINTS] += valid;
	}

	/* nothing more to produce: hold, and answer every request until told to drain (PROTOCOL.md §3.3) */
	ctx.state.store(HELD, std::memory_order_release);
	for (;;) {
		int st = ctx.state.load(std::memory_order_acquire);
		if (st == HOLD)
			ctx.state.store(HELD, std::memory_order_release);
		if (st == DRAIN) {
			ctx.state.store(QUIESCENT, std::memory_order_release);
			return;
		}
		cpu_relax();
	}
}

}
#endif
