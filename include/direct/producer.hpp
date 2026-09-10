#ifndef MITM_DIRECT_PRODUCER
#define MITM_DIRECT_PRODUCER

#include <algorithm>

#include "tools.hpp"
#include "router/router.hpp"
#include "direct/params.hpp"

/* The direct engine's producer: it evaluates the phase's function and pushes every image it gets. */

namespace mitm::direct {

/******************************** the producer thread ********************************/

/*
 * A producer's phase: evaluate the phase's function on its own contiguous piece of the phase's span -- the
 * round's chunk while filling, the whole domain while probing -- vlen inputs at a time, push every image as
 * (murmur64(image), preimage) to the receiver the hash's low word names, then close.  The pieces of the
 * n_producers producers cut the span into equal parts, the same on every node, so the partition needs no
 * message; the hash spares the producer a 64-bit division per point.
 */
template <class Wrapper>
void producer_round(Router_thread &rt, const Wrapper &wrapper, const Params &params, u64 *ctr, u64 round, int phase)
{
	constexpr int vlen = Wrapper::vlen;
	const u64 n_recv = Router_num_recv(rt);
	u64 lo = 0;
	u64 hi = params.domain;
	if (phase == FILL) {
		lo = round * params.per_round;
		hi = std::min(lo + params.per_round, params.domain);
	}
	u64 p = Router_rank(rt);
	u64 n_pieces = Router_num_send(rt);
	u64 span = hi - lo;
	u64 piece = span / n_pieces;
	u64 extra = span % n_pieces;
	u64 my_lo = lo + piece * p + std::min(p, extra);
	u64 my_hi = my_lo + piece + ((p < extra) ? 1 : 0);

	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	for (u64 base = my_lo; base < my_hi; base += vlen) {
		int valid = (my_hi - base < (u64) vlen) ? (int) (my_hi - base) : vlen;
		for (int k = 0; k < vlen; k++)
			x[k] = (k < valid) ? base + k : base;   /* a lane past the end still gets an input in the domain */
		wrapper.veval(phase, x, y);
		for (int k = 0; k < valid; k++) {
			u64 h = murmur64(y[k]);
			int dest = (int) (((h & 0xffffffffull) * n_recv) >> 32);
			Router_Push(h, x[k], dest, rt);
		}
		ctr[N_EVAL] += valid;
	}
	Router_Close(rt);
}

}
#endif
