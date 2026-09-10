#ifndef MITM_ROUTER_RING
#define MITM_ROUTER_RING

#include <atomic>
#include <cstddef>

#include "../tools.hpp"

namespace mitm {

/******************************** the ring ********************************/

/*
 * Lamport's single-producer / single-consumer ring of points: no
 * mutex, no CAS, each side owns one index and caches the other's.  Capacity rounds up to a
 * power of two.  Three cache lines: the consumer's, the producer's, and a read-only one for the geometry.
 * A receiver's inbox is one: block ids travel as points, `(block, count)`.
 */
class RouterRing {
public:
	/* public: Router_Dump reads head and tail from outside */
	alignas(64) Atomic<size_t> head;            /* written by the consumer only */
	size_t cached_tail;                        /* consumer-private copy of tail */
	alignas(64) Atomic<size_t> tail;            /* written by the producer only */
	size_t cached_head;                        /* producer-private copy of head */

private:
	alignas(64) const size_t capacity;         /* a power of two */
	const size_t mask;                         /* capacity - 1 */
	std::vector<Point> buf;                    /* the slots */

	static size_t round_up_pow2(size_t x)
	{
		size_t n = 1;
		while (n < x)
			n *= 2;
		return n;
	}

public:
	/* built by its owner: the zero-fill of `buf` is the first touch */
	RouterRing(size_t requested) : head(0), cached_tail(0), tail(0), cached_head(0),
	                               capacity(round_up_pow2(requested)), mask(capacity - 1), buf(capacity)
	{}

	/* producer side.  Returns false if the ring is full: the caller retries later, nothing is dropped. */
	bool push(const Point &x)
	{
		size_t t = tail;
		if (t - cached_head == capacity) {
			cached_head = head.load_acquire();
			if (t - cached_head == capacity)
				return false;                 /* really full */
		}
		buf[t & mask] = x;
		tail.store_release(t + 1);
		return true;
	}

	/* consumer side */
	bool pop(Point &x)
	{
		size_t h = head;
		if (h == cached_tail) {
			cached_tail = tail.load_acquire();
			if (h == cached_tail)
				return false;                 /* really empty */
		}
		x = buf[h & mask];
		head.store_release(h + 1);
		return true;
	}

	/* consumer side */
	bool empty()
	{
		size_t h = head;
		if (h == cached_tail)
			cached_tail = tail.load_acquire();
		return (h == cached_tail);
	}
};

}
#endif
