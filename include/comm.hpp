#ifndef MITM_COMM
#define MITM_COMM

#include <mpi.h>
#include <atomic>
#include <mutex>
#include <vector>
#include <array>
#include <memory>
#include <cmath>
#include <cassert>
#include <strings.h>          // ffsll, for the HyperLogLog

#include "tools.hpp"
#include "parameters.hpp"
#include "spsc.hpp"

namespace mitm {

class PcsDict;   /* inserter.hpp; complete in run(), where the one SharedContext is built and destroyed */


/****************************** collision queue *******************************/

struct CollisionCandidate {
	u64 i;                  /* mixing function version; a walker asserts it is the round's */
	u64 seed0;              /* the incoming point: its chain index, ... */
	u64 len0;               /* ... its trail length, ... */
	u64 end;                /* ... and its endpoint as the dictionary key, x / n_inserters */
	u64 seed1;              /* the point that was in the slot: its chain index, ... */
	u64 len1_maybe;         /* ... and its length; 0 == saturated in the dictionary, the walker re-walks it */
};

/*
 * Bounded, mutex-protected queue of dictionary hits: ONE per inserter, pushed by that inserter and
 * popped by the walkers of its thread group.  Candidates cross it in runs, never one at a time, so
 * that the lock and the `count` line are touched once per run.  Full == the run is truncated and the
 * inserter tallies the rest (DROP_COLL).  PROTOCOL.md §3.2.
 */
class CollisionQueue {
	std::mutex mtx;
	std::vector<CollisionCandidate> buf;
	size_t head = 0;                /* next to pop */
	size_t tail = 0;                /* next to push */
	std::atomic<size_t> count;      /* probed without the lock on the hot path */

public:
	CollisionQueue(size_t capacity) : buf(capacity < 1 ? 1 : capacity), count(0) {}

	/* Relaxed on purpose: a stale answer costs a walker one wasted pop() at worst. */
	bool is_empty() const
	{
		return count.load(std::memory_order_relaxed) == 0;
	}

	/* producer side.  Returns how many of the run fitted; the caller tallies the rest. */
	size_t push_bulk(const CollisionCandidate *in, size_t n)
	{
		std::lock_guard<std::mutex> lock(mtx);
		size_t room = buf.size() - count.load(std::memory_order_relaxed);
		size_t k = (room < n) ? room : n;
		for (size_t t = 0; t < k; t++) {
			buf[tail] = in[t];
			if (++tail == buf.size())
				tail = 0;
		}
		count.fetch_add(k, std::memory_order_release);
		return k;
	}

	/* consumer side.  Several walkers of the group may be here at once; the lock is what orders them. */
	size_t pop_bulk(CollisionCandidate *out, size_t max)
	{
		std::lock_guard<std::mutex> lock(mtx);
		size_t avail = count.load(std::memory_order_relaxed);
		size_t k = (avail < max) ? avail : max;
		for (size_t t = 0; t < k; t++) {
			out[t] = buf[head];
			if (++head == buf.size())
				head = 0;
		}
		count.fetch_sub(k, std::memory_order_release);
		return k;
	}
};


/********************************** counters **********************************/

/*
 * Every tally of the engine.  One u64[N_COUNTERS] per thread (ThreadContext::ctr) is also the wire
 * layout of a progress report (deltas) and of the end-of-round MPI_Reduce (totals): PROTOCOL.md §2.2.
 * None of this is part of the attack: it is what tells us whether the parameters are any good.
 */
enum counter {
	/* walkers */
	N_EVAL = 0,             /* evaluations of the mixing function, by walking or by resolving */
	N_DP,                   /* distinguished points found */
	N_POINTS_TRAILS,        /* sum of the lengths of the trails that reached a DP */
	N_COLLISIONS,           /* collisions located */
	COLLIDING_LEN_MIN,      /* sum of the shorter length of each colliding pair */
	COLLIDING_LEN_MAX,      /* ... and of the longer one */
	BAD_DP,                 /* trail gave up before a distinguished point, walking or re-walked to be measured */
	BAD_COLLISION,          /* the two trails "collide" on the same value */
	BAD_WALK_ROBINHOOD,     /* one trail is a suffix of the other */
	BAD_WALK_NONCOLLIDING,  /* dictionary false positive: the trails never meet */
	DROP_WALKERQ,           /* DP dropped: the queue to the comm thread was full */
	/* inserters */
	N_PROBE,                /* dictionary probes retired */
	BAD_PROBE,              /* dictionary slot was empty, or held a different key */
	DROP_COLL,              /* candidate dropped: the collision queue was full */
	/* comm thread */
	DROP_OUT,               /* DP dropped: the outgoing MPI buffer was still in flight */
	DROP_INSERTERQ,         /* DP dropped: a local inserter's queue was full */
	N_COUNTERS
};


/******************************** HyperLogLog *********************************/

static constexpr int HLL_REGISTERS = 0x10000;   /* one per value of the top 16 bits of a pair's hash */

/*
 * How many DISTINCT collisions a round found.  One HyperLogLog per walker, plain `u8`, written by its
 * owner alone and merged by the comm thread once every walker is quiescent: PROTOCOL.md §3.4.  None of
 * this is part of the attack -- it is what tells us whether the parameters are any good.
 */

/* record the collision (x0, x1): raise its register to the pair's rho, never lower it */
static inline void hll_record(u8 h[HLL_REGISTERS], u64 x0, u64 x1)
{
	u64 hash = murmur128(x0, x1);
	u64 idx = hash >> 48;
	u8 rho = (u8) ffsll(hash);
	if (h[idx] < rho)
		h[idx] = rho;
}

/* fold `src` into `dst`: the estimate below only ever needs the max of each register */
static inline void hll_merge(u8 dst[HLL_REGISTERS], const u8 src[HLL_REGISTERS])
{
	for (int k = 0; k < HLL_REGISTERS; k++)
		if (dst[k] < src[k])
			dst[k] = src[k];
}

/* the estimate, on merged registers (one round's, or the controller's all-time) */
static inline u64 distinct_collisions_estimation(const u8 h[HLL_REGISTERS])
{
	double acc = 0;
	double alpha = 0.7213 / (1 + 1.079 / HLL_REGISTERS);
	for (int i = 0; i < HLL_REGISTERS; i++)
		acc += std::ldexp(1.0, -(int) h[i]);       /* 2^-h[i]; (1 << h[i]) overflows past 31 */
	double E = alpha * ((double) HLL_REGISTERS * HLL_REGISTERS) / acc;
	if (E >= 2.5 * HLL_REGISTERS)
		return E;
	// low cardinality, potential correction
	int V = 0;
	for (int i = 0; i < HLL_REGISTERS; i++)
		if (h[i] == 0)
			V += 1;
	if (V == 0)
		return E;
	else
		return HLL_REGISTERS * std::log((double) HLL_REGISTERS / V);
}


/******************************** thread state ********************************/

/*
 * Winding a round down.  The comm thread writes a worker's state to ask, the worker writes it to
 * answer; the comm thread's own state is the phase of its round loop.  Values, transitions and the
 * order of the drain: PROTOCOL.md §3.3 and §4.5.
 */
enum thread_state {
	RUNNING, HOLD, HELD, DRAIN, QUIESCENT,                                /* workers (RUNNING, QUIESCENT: everyone) */
	COLLECTING, FLUSHING, WAITING, DRAINING_INSERTERS, DRAINING_WALKERS   /* comm thread only */
};


/******************************* thread context *******************************/

/*
 * What one thread owns that the comm thread also touches: its role, where it landed, its wind-down
 * state, its tallies and, for a worker, its SPSC queue.  One per thread, SharedContext::ctx[tid].
 */
struct alignas(64) ThreadContext {
	const int role;                 /* the comm thread dispatches on it */
	const int cpu;                  /* where the thread runs, asked to the kernel once pinned */
	const int numa_node;            /* its NUMA node: the one its first touch lands on */
	const int group;                /* its thread group; a walker's IS the inserter it resolves for (§1) */
	std::atomic<int> state;         /* PROTOCOL.md §3.3.  Thread 0's is the phase of its round loop */
	u64 ctr[N_COUNTERS] = {};       /* this thread's alone, plain u64; exact only once it is QUIESCENT (§3.4) */
	std::unique_ptr<SPSCQueue> q;   /* walker: to the comm thread; inserter: from it; null for the comm thread */
	std::vector<u8> hll;            /* walker: its own HyperLogLog of the round's collisions; empty otherwise */

	/* built by the thread it describes, once pinned: its buffers are then its first touch (§4.1) */
	ThreadContext(int role, int cpu, int numa_node, int group, const Parameters &params)
		: role(role), cpu(cpu), numa_node(numa_node), group(group), state(RUNNING)
	{
		if (role == WALKER) {
			q = std::make_unique<SPSCQueue>(params.walker_queue_capacity);
			hll.assign(HLL_REGISTERS, 0);
		} else if (role == INSERTER)
			q = std::make_unique<SPSCQueue>(params.inserter_queue_capacity);
	}
};


/******************************* shared context *******************************/

/*
 * What every thread of the node shares.  One per rank, a local of run(); ctx, shards and coll_q are
 * sized here and filled by each thread once it is pinned (PROTOCOL.md §4.1).
 */
struct SharedContext {
	/* the round header: written by the comm thread, read by everyone after an omp barrier */
	u64 i = 0;                      /* mixing function version */
	u64 root_seed = 0;              /* chain j starts at root_seed + j * multiplier */
	u64 stop = 0;                   /* no next round: every thread leaves the loop */

	/* one per thread, by tid; a worker only ever sees its own */
	std::vector<std::unique_ptr<ThreadContext>> ctx;

	/* shard r: inserter r (thread 1 + r) alone builds, probes and flushes it */
	std::vector<std::unique_ptr<PcsDict>> shards;

	/* queue r: inserter r pushes, the walkers of group r pop (PROTOCOL.md §3.2) */
	std::vector<std::unique_ptr<CollisionQueue>> coll_q;

	/* the golden pair: any walker writes (the first one wins), the comm thread polls */
	std::mutex golden_mtx;
	std::atomic<bool> golden_found;
	u64 golden[3];              /* i, x0, x1: the TAG_SOLUTION payload (solution_field) */

	SharedContext(const Parameters &params)
		: ctx(params.n_threads), shards(params.inserters_per_node),
		  coll_q(params.inserters_per_node), golden_found(false)
	{
		golden[0] = golden[1] = golden[2] = 0;
	}

	void set_golden(u64 i, u64 x0, u64 x1)
	{
		std::lock_guard<std::mutex> lock(golden_mtx);
		if (golden_found.load(std::memory_order_relaxed))
			return;                        /* keep the first one */
		golden[0] = i;
		golden[1] = x0;
		golden[2] = x1;
		golden_found.store(true, std::memory_order_release);
	}
};


/*************************** outgoing bulk DP buffers *************************/

/*
 * The outgoing bulk DP messages: one double buffer per destination, `ready` filling while `outgoing`
 * is in flight.  Never blocks: a point that finds both busy is dropped and push() says so (the caller
 * tallies DROP_OUT).  A send is tested only when its slot is wanted: PROTOCOL.md §2.1.
 */
class OutBuffers {
	using Buffer = std::vector<u64>;

	MPI_Comm comm;
	int n_nodes;
	size_t cap;                              /* u64 per buffer */
	std::vector<Buffer> ready;               /* being filled */
	std::vector<Buffer> outgoing;            /* in flight */
	std::vector<MPI_Request> req;            /* for the OUTGOING buffers */

	/* move `ready` out, if the previous send is done.  false == still busy */
	bool rotate(int dst)
	{
		if (req[dst] != MPI_REQUEST_NULL) {
			int flag = 0;
			MPI_Test(&req[dst], &flag, MPI_STATUS_IGNORE);
			if (!flag)
				return false;
		}
		outgoing[dst].clear();
		std::swap(ready[dst], outgoing[dst]);
		if (not outgoing[dst].empty())       /* never send empty: that is the sentinel */
			MPI_Isend(outgoing[dst].data(), outgoing[dst].size(), MPI_UINT64_T, dst, TAG_POINTS, comm, &req[dst]);
		return true;
	}

public:
	OutBuffers(MPI_Comm comm, int n_nodes, size_t dp_capacity)
		: comm(comm), n_nodes(n_nodes), cap(DP_WORDS * dp_capacity),
		  ready(n_nodes), outgoing(n_nodes), req(n_nodes, MPI_REQUEST_NULL)
	{
		for (int d = 0; d < n_nodes; d++) {
			ready[d].reserve(cap);
			outgoing[d].reserve(cap);
		}
	}

	/* append one DP to the buffer bound for node `dst`.  false == dropped. */
	bool push(const DP &p, int dst)
	{
		if ((ready[dst].size() + DP_WORDS > cap) && (not rotate(dst)))
			return false;
		ready[dst].push_back(p.seed);
		ready[dst].push_back(p.x);
		ready[dst].push_back(p.len);
		return true;
	}

	/*
	 * End of round: push out what is left and test the sends.  True once nothing is left and nothing
	 * is in flight.  Never waits: call it from a loop that keeps receiving, or two nodes deadlock.
	 */
	bool flush_poll()
	{
		bool done = true;
		for (int d = 0; d < n_nodes; d++) {
			if (not ready[d].empty()) {
				rotate(d);
				done = false;
			} else if (req[d] != MPI_REQUEST_NULL) {
				int flag = 0;
				MPI_Test(&req[d], &flag, MPI_STATUS_IGNORE);
				if (!flag)
					done = false;
			}
		}
		return done;
	}

	/*
	 * Tell every node, self included, that our round's points are all sent: one zero-length TAG_POINTS
	 * each, by MPI_Bsend (PROTOCOL.md §2.1).  Needs the engine's Bsend buffer, attached by run() first.
	 */
	void send_sentinels()
	{
		for (int d = 0; d < n_nodes; d++)
			MPI_Bsend(NULL, 0, MPI_UINT64_T, d, TAG_POINTS, comm);
	}
};

/* the TAG_SOLUTION payload and SharedContext::golden, word by word (PROTOCOL.md §2.2) */
enum solution_field {
	SOL_I = 0, SOL_X0, SOL_X1,
	SOL_NWORDS
};

}
#endif
