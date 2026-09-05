#ifndef MITM_COMM
#define MITM_COMM

#include <mpi.h>
#include <atomic>
#include <mutex>
#include <vector>
#include <memory>
#include <type_traits>
#include <cassert>

#include "tools.hpp"
#include "parameters.hpp"
#include "spsc.hpp"

namespace mitm {

/******************************** thread state ********************************/

/*
 * Winding a round down.  The comm thread writes a worker's state to ask, the worker writes it to
 * answer; the comm thread's own state is the phase of its round loop.  Values, transitions and the
 * order of the drain: PROTOCOL.md §3.3 and §4.5.
 */
enum thread_state {
	RUNNING, HOLD, HELD, DRAIN, QUIESCENT,                                /* workers (RUNNING, QUIESCENT: everyone) */
	COLLECTING, FLUSHING, WAITING, DRAINING_DICTS, DRAINING_PRODUCERS     /* comm thread only */
};


/******************************* thread context *******************************/

/*
 * What one thread owns that the comm thread also touches: its role, where it landed, its wind-down
 * state, its tallies, for a worker its SPSC queue, and the scheme's own per-thread statistics (PCS: a
 * producer's HyperLogLog).  One per thread, SharedContext::ctx[tid].
 */
template <class Scheme>
struct alignas(64) ThreadContext {
	const int role;                 /* the comm thread dispatches on it */
	const int cpu;                  /* where the thread runs, asked to the kernel once pinned */
	const int numa_node;            /* its NUMA node: the one its first touch lands on */
	const int group;                /* its thread group; a producer's IS the dict thread it works for (§1) */
	std::atomic<int> state;         /* PROTOCOL.md §3.3.  Thread 0's is the phase of its round loop */
	u64 ctr[Scheme::N_COUNTERS] = {};   /* this thread's alone, plain u64; exact once it is QUIESCENT (§3.4) */
	std::unique_ptr<SPSCQueue> q;   /* producer: to the comm thread; dict thread: from it; null for the comm thread */
	typename Scheme::ThreadStats scheme;   /* the scheme's round-scoped statistics of this thread (§3.4) */

	/* built by the thread it describes, once pinned: its buffers are then its first touch (§4.1) */
	ThreadContext(int role, int cpu, int numa_node, int group, const Parameters &params)
		: role(role), cpu(cpu), numa_node(numa_node), group(group), state(RUNNING), scheme(role)
	{
		if (role == PRODUCER)
			q = std::make_unique<SPSCQueue>(params.producer_queue_capacity);
		else if (role == DICT)
			q = std::make_unique<SPSCQueue>(params.dict_queue_capacity);
	}
};


/******************************* shared context *******************************/

/*
 * What every thread of the node shares.  One per rank, a local of run(); ctx and shards are sized here
 * and filled by each thread once it is pinned (PROTOCOL.md §4.1).  `scheme` is the scheme's own shared
 * state (PCS: its collision queues).
 */
template <class Scheme>
struct SharedContext {
	/* the round header: written by the comm thread, read by everyone after an omp barrier (§4.2) */
	typename Scheme::Header header;
	u64 stop = 0;                   /* no next round: every thread leaves the loop */

	/* one per thread, by tid; a worker only ever sees its own */
	std::vector<std::unique_ptr<ThreadContext<Scheme>>> ctx;

	/* shard r: dict thread r (thread 1 + r) alone builds, probes and flushes it */
	std::vector<std::unique_ptr<typename Scheme::Dict>> shards;

	/* what the scheme adds; built here, filled by the threads it belongs to */
	typename Scheme::Shared scheme;

	/* the golden pair: any worker writes (the first one wins), the comm thread polls */
	std::mutex golden_mtx;
	std::atomic<bool> golden_found;
	u64 golden[3];              /* i, x0, x1: the TAG_SOLUTION payload (solution_field) */

	SharedContext(const typename Scheme::Params &params)
		: ctx(params.n_threads), shards(params.dicts_per_node), scheme(params), golden_found(false)
	{
		static_assert(std::is_trivially_copyable<typename Scheme::Header>::value
		              && sizeof(typename Scheme::Header) % sizeof(u64) == 0,
		              "a round header goes on the wire as whole u64 words");
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


/************************** outgoing bulk point buffers ***********************/

/*
 * The outgoing bulk point messages: one double buffer per destination, `ready` filling while
 * `outgoing` is in flight.  Never blocks: a point that finds both busy is refused and push() says so
 * (the caller drops or holds it, PROTOCOL.md §5).  A send is tested only when its slot is wanted (§2.1).
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
	OutBuffers(MPI_Comm comm, int n_nodes, size_t point_capacity)
		: comm(comm), n_nodes(n_nodes), cap(POINT_WORDS * point_capacity),
		  ready(n_nodes), outgoing(n_nodes), req(n_nodes, MPI_REQUEST_NULL)
	{
		for (int d = 0; d < n_nodes; d++) {
			ready[d].reserve(cap);
			outgoing[d].reserve(cap);
		}
	}

	/* append one point to the buffer bound for node `dst`.  false == refused, nothing was written. */
	bool push(const Point &p, int dst)
	{
		if ((ready[dst].size() + POINT_WORDS > cap) && (not rotate(dst)))
			return false;
		ready[dst].push_back(p.key);
		ready[dst].push_back(p.val);
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
