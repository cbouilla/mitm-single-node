#ifndef MITM_SPSC
#define MITM_SPSC

#include <atomic>
#include <vector>
#include <cstddef>

#include "parameters.hpp"

namespace mitm {

/*
 * Single-producer / single-consumer bounded queue of distinguished points
 * (Lamport 1983), wait-free on both sides.  No mutex, no CAS: the producer owns
 * `tail`, the consumer owns `head`, and each keeps a private cached copy of the
 * other's index so that the common path never reads the other side's cache line.
 *
 * Capacity is rounded up to a power of two so that wrapping is a mask.
 *
 * One instance sits between each walker thread and the comm thread, and between
 * the comm thread and each inserter thread.
 */
class SPSCQueue {
public:
	/* each index sits alone in its own cache line.  Public so a third party (the
	   comm thread, at end of round) can test head == tail without an accessor. */
	alignas(64) std::atomic<size_t> head;      /* written by the consumer only */
	alignas(64) std::atomic<size_t> tail;      /* written by the producer only */

private:
	size_t capacity;
	size_t mask;
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
	/* init order follows declaration order: head/tail, then capacity (mask and buf
	   depend on it), then the private caches */
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
