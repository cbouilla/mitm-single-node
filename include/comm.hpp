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

class PcsDict;   /* inserter.hpp.  SharedContext::shards holds pointers to the shards, and the one
                    SharedContext is built and destroyed in run(), where the class is complete. */

/*
 * Where state lives, per rank.  Private to a worker thread and never touched by the
 * comm thread: a local of that thread's function.  Specific to one thread but also
 * touched by the comm thread: its ThreadContext.  Common to all threads: the
 * SharedContext.  Private to the comm thread: a member or a local of CommThread
 * (engine.hpp).  Two refinements: every counter, the comm thread's own included,
 * is in a ThreadContext, so that the node's tallies are one loop with no special
 * case; and the dictionary shards are in the SharedContext although only their
 * inserter touches them today, because zeroing them may one day be collective.
 */


/*
 * Routing.  The shard of a DP is `x % n_inserters`, with
 * n_inserters == n_nodes * inserters_per_node.  Because n_nodes divides n_inserters
 * that splits cleanly, and each piece is computed inline at the stage that needs it,
 * which keeps two of the three divisions off the comm thread (the only
 * single-threaded stage):
 *
 *     node = x % n_nodes                         (sending comm thread)
 *     slot = (x / n_nodes) % inserters_per_node  (receiving comm thread)
 *     key  = x / n_inserters                     (inserter thread)
 */


/****************************** collision queue *******************************/

struct CollisionCandidate {
	u64 i;                  /* mixing function index; stale entries are discarded */
	u64 seed0, len0, end;
	u64 seed1, len1_maybe;  /* 0 == length unknown --> walk_nolen1 */
};

/*
 * Bounded, mutex-protected.  Inserter threads push; walker threads pop and pay for
 * the walk.  A full queue drops the candidate and tells the pusher, which keeps the
 * tally (DROP_COLL in its ThreadContext) -- the queue itself counts nothing.  The
 * emptiness probe is a relaxed atomic so walkers stay off the lock on the common
 * path.
 */
class CollisionQueue {
	std::mutex mtx;
	std::vector<CollisionCandidate> buf;
	size_t head = 0, tail = 0;
	std::atomic<size_t> count;      /* probed without the lock on the hot path */

public:
	CollisionQueue(size_t capacity) : buf(capacity < 1 ? 1 : capacity), count(0) {}

	/* Relaxed on purpose: a stale answer costs a walker one wasted pop() at worst. */
	bool is_empty() const
	{
		return count.load(std::memory_order_relaxed) == 0;
	}

	bool push(const CollisionCandidate &c)
	{
		std::lock_guard<std::mutex> lock(mtx);
		if (count.load(std::memory_order_relaxed) == buf.size())
			return false;                /* full: the caller tallies the drop */
		buf[tail] = c;
		if (++tail == buf.size())
			tail = 0;
		count.fetch_add(1, std::memory_order_release);
		return true;
	}

	bool pop(CollisionCandidate &c)
	{
		std::lock_guard<std::mutex> lock(mtx);
		if (count.load(std::memory_order_relaxed) == 0)
			return false;
		c = buf[head];
		if (++head == buf.size())
			head = 0;
		count.fetch_sub(1, std::memory_order_release);
		return true;
	}
};


/********************************** counters **********************************/

/*
 * Every diagnostic tally of the engine, by index.  One array of u64 per thread
 * (ThreadContext::ctr) holds them all, whichever thread produces them, so that the
 * same layout serves for the per-thread tallies, the progress reports (deltas,
 * TAG_REPORT) and the end-of-round MPI_Reduce.
 *
 * None of this is part of the attack: it is what tells us whether the parameters are
 * any good.
 */
enum counter {
	/* walkers */
	N_EVAL = 0,             /* evaluations of the mixing function, by walking or by resolving */
	N_DP,                   /* distinguished points found */
	N_POINTS_TRAILS,        /* sum of the lengths of the trails that reached a DP */
	N_COLLISIONS,           /* collisions located */
	COLLIDING_LEN_MIN,      /* sum of the shorter length of each colliding pair */
	COLLIDING_LEN_MAX,      /* ... and of the longer one */
	BAD_DP,                 /* trail gave up before reaching a distinguished point */
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


/******************************** thread state ********************************/

/*
 * How a round is wound down.  Each worker thread carries its own state; the comm
 * thread drives the sequence and reads back the states the workers set for
 * themselves, so a thread always announces its own completion and there is nothing
 * to race on.
 *
 *   walker:   RUNNING -[comm]-> HOLD -[self]-> HELD -[comm]-> DRAIN -[self]-> QUIESCENT
 *   inserter: RUNNING ---------------------------[comm]-> DRAIN -[self]-> QUIESCENT
 *
 * HOLD tells a walker to stop producing points; HELD is its acknowledgement, which
 * the comm thread needs because a walker only looks at its state once per chunk --
 * an "empty" walker queue mid-chunk is not an idle one.  DRAIN means "nothing more
 * will ever arrive for you: finish what you hold and go quiet".
 *
 * The comm thread's own state is the phase of its round loop (CommThread::comm_round
 * in engine.hpp), which it alone reads and writes:
 *
 *   comm:     RUNNING -> COLLECTING -> FLUSHING -> WAITING -> DRAINING_INSERTERS
 *                     -> DRAINING_WALKERS -> QUIESCENT
 *
 * The order is what makes a round airtight: no thread declares itself finished while
 * something can still arrive for it.  RUNNING is the steady state, until the
 * end-of-round signal; then the walkers are told to HOLD.  COLLECTING waits for them
 * all to be HELD and routes what they left in their queues.  FLUSHING ships the
 * partial output buffers and, once every send has completed, tells every node
 * (itself included) that we are done sending.  WAITING takes delivery until every
 * node has said the same, so nothing can arrive for an inserter any more, and the
 * inserters are told to DRAIN.  DRAINING_INSERTERS waits for them to be QUIESCENT:
 * no new collision candidate can appear, so the walkers are told to DRAIN in turn.
 * DRAINING_WALKERS waits for the walkers to be QUIESCENT -- last, which is what keeps
 * a golden pair found on the very last candidate from being lost.  Then the comm
 * thread is QUIESCENT itself and its round is over.
 */
enum thread_state {
	RUNNING, HOLD, HELD, DRAIN, QUIESCENT,                                /* workers (RUNNING, QUIESCENT: everyone) */
	COLLECTING, FLUSHING, WAITING, DRAINING_INSERTERS, DRAINING_WALKERS   /* comm thread only */
};


/******************************* thread context *******************************/

/*
 * What one thread owns that the comm thread also touches: its role, its wind-down
 * state, its tallies and, for a worker, its SPSC queue.  One per thread, the comm
 * thread's included (it has no queue; its tallies are here like everyone's, and its
 * state is the phase of its round loop).
 */
struct alignas(64) ThreadContext {
	const int role;                 /* the comm thread dispatches on it */
	const int cpu;                  /* where the thread runs, asked to the kernel once pinned */
	const int numa_node;            /* its NUMA node: the one its first touch lands on */

	/* see thread_state above: the comm thread asks, the thread answers.  Thread 0's
	   holds the phase of its own round loop, which nobody else reads. */
	std::atomic<int> state;

	/* Every tally of this thread, written by it alone with no synchronisation
	   (plain u64, see enum counter).  The comm thread sums them over the threads
	   after an `omp flush` for its progress reports, where a little staleness does
	   not matter, and reads them exactly once every worker is QUIESCENT (the store
	   that says so is a release, and it loads it with acquire), when it also clears
	   them. */
	u64 ctr[N_COUNTERS] = {};

	/* the thread's own queue.  Null for the comm thread, which has none. */
	std::unique_ptr<SPSCQueue> q;      /* walker: to the comm thread.  inserter: from it */

	/* Built by the thread it describes, once that thread is pinned (see run() in
	   engine.hpp): the queue's buffer is then that thread's first touch, which is
	   what places it on its NUMA node.  `cpu` and `numa_node` are what the kernel
	   answered after pinning; with --no-bind, only where the thread happened to be.
	   Shard r's NUMA node is ctx[1 + r]->numa_node. */
	ThreadContext(int role, int cpu, int numa_node, const Parameters &params)
		: role(role), cpu(cpu), numa_node(numa_node), state(RUNNING)
	{
		if (role == WALKER)
			q = std::make_unique<SPSCQueue>(params.walker_queue_capacity);
		else if (role == INSERTER)
			q = std::make_unique<SPSCQueue>(params.inserter_queue_capacity);
	}
};


/******************************* shared context *******************************/

static constexpr int HLL_REGISTERS = 0x10000;

/*
 * What every thread of the node shares: the round header, the per-thread table, the
 * dictionary shards, the collision queue, the HyperLogLog and the golden pair.  One
 * per rank, a local of run(), built before the team.  The two tables are sized here
 * and filled inside the team by each thread once it is pinned (NUMA first touch, see
 * run()); the shards are here rather than in their inserter's ThreadContext because
 * zeroing them may one day be collective.
 */
struct SharedContext {
	/* round header: written by the comm thread, read by everyone after an omp barrier */
	u64 i = 0;
	u64 root_seed = 0;
	u64 stop = 0;

	/* one slot per thread, indexed by tid.  Only the comm thread and run() index it;
	   a worker is handed its own slot by reference and never looks at the others. */
	std::vector<std::unique_ptr<ThreadContext>> ctx;

	/* the dictionary shards: shard r belongs to inserter r, thread 1 + r, which alone
	   builds, probes and flushes it (today) */
	std::vector<std::unique_ptr<PcsDict>> shards;

	/* inserters push, walkers pop */
	alignas(64) CollisionQueue coll_q;

	/*
	 * HyperLogLog over this round's collisions, one per node, used to estimate how
	 * many of them are DISTINCT -- the quantity that actually drives the attack.  Any
	 * walker raises a register with a compare-and-swap (an atomic max); nobody reads
	 * them before every walker is QUIESCENT, when the comm thread copies them out for
	 * the MPI_MAX reduction and zeroes them.  Relaxed throughout, for that reason: the
	 * walkers' QUIESCENT stores (release) order them for the comm thread's acquire.
	 * C++17: a default-constructed std::atomic is NOT zero, hence the loop in the
	 * constructor.
	 */
	alignas(64) std::atomic<u8> hll[HLL_REGISTERS];
	static_assert(std::atomic<u8>::is_always_lock_free, "the HyperLogLog registers must be lock-free");

	/* the golden pair: any walker writes (the first one wins), the comm thread polls */
	std::mutex golden_mtx;
	std::atomic<bool> golden_found;
	u64 golden[3];              /* i, x0, x1: the TAG_SOLUTION payload (solution_field) */

	SharedContext(const Parameters &params)
		: ctx(params.n_threads), shards(params.inserters_per_node),
		  coll_q(params.coll_queue_capacity), golden_found(false)
	{
		for (int k = 0; k < HLL_REGISTERS; k++)
			hll[k].store(0, std::memory_order_relaxed);
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

	/*
	 * A collision (x0, x1), x0 < x1: raise its HyperLogLog register.  The register is
	 * the top 16 bits of the pair's hash, its value the position of the lowest set
	 * bit.  The CAS loop only ever raises a register, and stops as soon as it sees
	 * one at least as high.
	 */
	void found_collision(u64 x0, u64 x1)
	{
		u64 h = murmur128(x0, x1);
		u64 idx = h >> 48;
		u8 rho = (u8) ffsll(h);
		u8 cur = hll[idx].load(std::memory_order_relaxed);
		while (cur < rho && not hll[idx].compare_exchange_weak(cur, rho, std::memory_order_relaxed))
			;                              /* a failure reloaded cur: go round again */
	}

	/* the HyperLogLog estimator, on a plain copy of the registers (one round's, or the
	   controller's all-time) */
	static u64 distinct_collisions_estimation(const u8 h[HLL_REGISTERS])
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
};


/*************************** outgoing bulk DP buffers *************************/

/*
 * Double buffering, one pair per destination node: `ready` accumulates, `outgoing`
 * is in flight.  Never blocks -- if the previous send has not finished when `ready`
 * fills up, the point is dropped, which costs work but never correctness; push()
 * says so and the comm thread keeps the tally (DROP_OUT), the buffers count nothing.
 *
 * Because the transmitting data lives in `outgoing` and push() only ever touches
 * `ready`, a buffer that is still in flight simply cannot be handed to MPI twice.
 *
 * The sends are never polled for their own sake: rotate() tests a request exactly
 * when its slot is wanted, and flush_poll() at the end of the round.  MPI progress
 * is per process, not per request -- every MPI_Test the comm loop makes on its
 * receives drives these sends too -- so a periodic MPI_Testsome over n_nodes
 * requests would only cost a pass over them on the single-threaded stage.
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
				return false;                /* the outgoing buffer is still being sent */
		}
		/* The outgoing buffer is available */
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
	 * End of round: get every remaining point out.  Returns true once nothing is left
	 * to send and nothing is still in flight.  Call it in a loop that also drains
	 * incoming traffic -- blocking here deadlocks, because two nodes both sitting in
	 * MPI_Wait stop reposting their receives and neither one's sends can land.
	 */
	bool flush_poll()
	{
		bool done = true;
		for (int d = 0; d < n_nodes; d++) {
			if (not ready[d].empty()) {
				rotate(d);                   /* may have to wait for the previous send */
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
	 * Zero-length message == "I have sent you everything for this round".
	 * Ordering holds -- MPI messages between a pair of ranks are
	 * non-overtaking whatever the send mode, so a sentinel initiated after
	 * the data Isends is delivered after them.
	 *
	 * Requires the engine's MPI_Bsend buffer to be attached; run() does that before
	 * anything else of the engine exists, long before any round drains.
	 */
	void send_sentinels()
	{
		for (int d = 0; d < n_nodes; d++)
			MPI_Bsend(NULL, 0, MPI_UINT64_T, d, TAG_POINTS, comm);
	}
};

enum solution_field {
	SOL_I = 0, SOL_X0, SOL_X1,
	SOL_NWORDS
};

}
#endif
