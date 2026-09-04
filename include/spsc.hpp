#ifndef MITM_SPSC
#define MITM_SPSC

#include <atomic>
#include <vector>
#include <cstddef>

namespace mitm {

/*
 * Lamport's single-producer / single-consumer ring of DPs: no mutex, no CAS, each side owns one index
 * and caches the other's (PROTOCOL.md §3.1).  Capacity rounds up to a power of two.  One between each
 * walker and the comm thread, one between the comm thread and each inserter.
 */
class SPSCQueue {
public:
	/* public: the comm thread tests head == tail at the end of a round, without an accessor */
	alignas(64) std::atomic<size_t> head;      /* written by the consumer only */
	alignas(64) std::atomic<size_t> tail;      /* written by the producer only */

private:
	size_t capacity;                           /* a power of two */
	size_t mask;                               /* capacity - 1 */
	std::vector<DP> buf;

	alignas(64) size_t cached_head;            /* producer-private copy of head */
	alignas(64) size_t cached_tail;            /* consumer-private copy of tail */

	static size_t round_up_pow2(size_t x)
	{
		size_t n = 1;
		while (n < x)
			n *= 2;
		return n;
	}

public:
	SPSCQueue(size_t requested) : head(0), tail(0), capacity(round_up_pow2(requested)),
	                              mask(capacity - 1), buf(capacity),
	                              cached_head(0), cached_tail(0)
	{}

	/* producer side.  Returns false if the queue is full (the caller drops the item). */
	bool push(const DP &x)
	{
		size_t t = tail.load(std::memory_order_relaxed);
		if (t - cached_head == capacity) {
			cached_head = head.load(std::memory_order_acquire);
			if (t - cached_head == capacity)
				return false;                 /* really full */
		}
		buf[t & mask] = x;
		tail.store(t + 1, std::memory_order_release);
		return true;
	}

	/* consumer side */
	bool pop(DP &x)
	{
		size_t h = head.load(std::memory_order_relaxed);
		if (h == cached_tail) {
			cached_tail = tail.load(std::memory_order_acquire);
			if (h == cached_tail)
				return false;                 /* really empty */
		}
		x = buf[h & mask];
		head.store(h + 1, std::memory_order_release);
		return true;
	}

	/* consumer side.  Grab up to `max` items at once, amortizing the atomics. */
	size_t pop_bulk(DP *out, size_t max)
	{
		size_t h = head.load(std::memory_order_relaxed);
		if (h == cached_tail) {
			cached_tail = tail.load(std::memory_order_acquire);
			if (h == cached_tail)
				return 0;
		}
		size_t avail = cached_tail - h;
		size_t n = (avail < max) ? avail : max;
		for (size_t k = 0; k < n; k++)
			out[k] = buf[(h + k) & mask];
		head.store(h + n, std::memory_order_release);
		return n;
	}

	/* consumer side */
	bool empty()
	{
		size_t h = head.load(std::memory_order_relaxed);
		if (h == cached_tail)
			cached_tail = tail.load(std::memory_order_acquire);
		return (h == cached_tail);
	}


};

}
#endif
