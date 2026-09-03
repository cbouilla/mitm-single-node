#ifndef MITM_COMM
#define MITM_COMM

#include <mpi.h>
#include <atomic>
#include <mutex>
#include <vector>
#include <array>
#include <memory>

#include "counters.hpp"
#include "parameters.hpp"
#include "spsc.hpp"
#include "dict.hpp"

namespace mitm {


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
 * tally (DROP_COLL in its Counters) -- the queue itself counts nothing.  The
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


/******************************** round state *********************************/

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
 */
enum thread_state {RUNNING, HOLD, HELD, DRAIN, QUIESCENT};

struct RoundState {
	/* written by the comm thread, read by everyone after an omp barrier */
	u64 i = 0;
	u64 root_seed = 0;
	u64 stop = 0;


	std::mutex golden_mtx;
	std::atomic<bool> golden_found;
	u64 golden[3];              /* i, x0, x1: the TAG_SOLUTION payload (solution_field) */

	RoundState() : golden_found(false)
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


/******************************* thread context *******************************/

/*
 * What a worker thread owns -- its SPSC queue and, for an inserter, its dictionary
 * shard -- plus the few fields the comm thread also touches.  One per thread, the
 * comm thread's included (it owns neither a queue nor a shard: its objects are the
 * CommThread's; what it does own here is its tallies).
 */
struct alignas(64) ThreadContext {
	const int role;                 /* the comm thread dispatches on it */
	const int cpu;                  /* where the thread pinned itself; rank 0 prints the map */

	/* see thread_state above: the comm thread asks, the thread answers */
	std::atomic<int> state;

	/* Every tally of this thread, written by it alone with no synchronisation
	   (plain u64).  The comm thread sums them over the threads after an `omp flush`
	   for its progress reports, where a stale unit or two does not matter, and reads
	   them exactly after the end-of-round barrier, when it also clears them. */
	Counters ctr;

	/* the thread's own.  Null for the roles that have none: the comm thread has no
	   queue, only an inserter has a shard.  Read-only pointers once built. */
	std::unique_ptr<SPSCQueue> q;      /* walker: to the comm thread.  inserter: from it */
	std::unique_ptr<PcsDict> dict;     /* inserter: its shard of the dictionary */

	/* Built by the thread it describes, once that thread is pinned (see run() in engine.hpp):
	   the queue's buffer and the shard's zero-fill are then that thread's first touch,
	   which is what places them on its NUMA node. */
	ThreadContext(int role, int cpu, const Parameters &params)
		: role(role), cpu(cpu), state(RUNNING)
	{
		if (role == WALKER) {
			q = std::make_unique<SPSCQueue>(params.walker_queue_capacity);
		} else if (role == INSERTER) {
			q = std::make_unique<SPSCQueue>(params.inserter_queue_capacity);
			dict = std::make_unique<PcsDict>(params.jbits, params.w_shard);
		}
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
 */
class OutBuffers {
	using Buffer = std::vector<u64>;

	MPI_Comm comm;
	int n_nodes;
	size_t cap;                              /* u64 per buffer */
	std::vector<Buffer> ready;               /* being filled */
	std::vector<Buffer> outgoing;            /* in flight */
	std::vector<MPI_Request> req;            /* for the OUTGOING buffers */
	std::vector<int> done_idx;

	void start_send(int dst)
	{
		if (outgoing[dst].empty())
			return;                          /* never send empty: that is the sentinel */
		MPI_Isend(outgoing[dst].data(), outgoing[dst].size(), MPI_UINT64_T, dst, TAG_POINTS,
		          comm, &req[dst]);
		bytes_sent += outgoing[dst].size() * sizeof(u64);
	}

	/* move `ready` out, if the previous send is done.  false == still busy */
	bool rotate(int dst)
	{
		if (req[dst] != MPI_REQUEST_NULL) {
			int flag = 0;
			MPI_Test(&req[dst], &flag, MPI_STATUS_IGNORE);
			if (!flag)
				return false;
		}
		outgoing[dst].clear();               /* clear BEFORE the swap: leaves ready empty */
		std::swap(ready[dst], outgoing[dst]);
		start_send(dst);
		return true;
	}

public:
	u64 bytes_sent = 0;

	OutBuffers(MPI_Comm comm, int n_nodes, size_t dp_capacity)
		: comm(comm), n_nodes(n_nodes), cap(DP_WORDS * dp_capacity),
		  ready(n_nodes), outgoing(n_nodes), req(n_nodes, MPI_REQUEST_NULL), done_idx(n_nodes)
	{
		for (int d = 0; d < n_nodes; d++) {
			ready[d].reserve(cap);
			outgoing[d].reserve(cap);
		}
	}

	/* append one DP to the buffer bound for node `dst`.  false == dropped. */
	bool push(const DP &p, int dst)
	{
		if (ready[dst].size() + DP_WORDS > cap && not rotate(dst))
			return false;                    /* the caller tallies the drop */
		ready[dst].push_back(p.seed);
		ready[dst].push_back(p.x);
		ready[dst].push_back(p.len);
		return true;
	}

	/* release whatever finished transmitting */
	void poll()
	{
		int outcount = 0;
		MPI_Testsome(req.size(), req.data(), &outcount, done_idx.data(), MPI_STATUSES_IGNORE);
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
	 * Zero-length message == "I have sent you everything for this round".  Buffered
	 * mode: it completes locally, so there is no request to track and nothing to wait
	 * on, and a zero-length message costs only MPI_BSEND_OVERHEAD of the attached
	 * buffer.  Ordering still holds -- MPI messages between a pair of ranks are
	 * non-overtaking whatever the send mode, so a sentinel initiated after the data
	 * Isends is delivered after them.
	 *
	 * Requires the ControlChannel's MPI_Buffer_attach to have happened already; it
	 * is done at construction, long before any round drains.
	 */
	void send_sentinels()
	{
		for (int d = 0; d < n_nodes; d++)
			MPI_Bsend(NULL, 0, MPI_UINT64_T, d, TAG_POINTS, comm);
	}
};


/*************************** incoming bulk DP buffers *************************/

/*
 * A pool of MPI_ANY_SOURCE receives, so buffer memory is set by `n_in_buffers`
 * rather than by the number of peers.
 *
 * Deliberately just state: the comm thread drives the Testsome/scatter/repost loop
 * itself (CommThread::poll_incoming), because handing this class a scatter callback
 * would mean a closure at every call site.
 */
struct InBuffers {
	MPI_Comm comm;
	size_t cap;                                 /* u64 per buffer */
	std::vector<std::vector<u64>> data;
	std::vector<MPI_Request> req;
	std::vector<int> done_idx;
	std::vector<MPI_Status> done_st;

	int n_sentinels = 0;
	u64 bytes_recv = 0;

	InBuffers(MPI_Comm comm, int n_buffers, size_t dp_capacity)
		: comm(comm), cap(DP_WORDS * dp_capacity), data(n_buffers),
		  req(n_buffers, MPI_REQUEST_NULL), done_idx(n_buffers), done_st(n_buffers)
	{
		for (int k = 0; k < n_buffers; k++) {
			data[k].resize(cap);
			MPI_Irecv(data[k].data(), cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS, comm, &req[k]);
		}
	}

	void shutdown()
	{
		for (size_t k = 0; k < req.size(); k++)
			if (req[k] != MPI_REQUEST_NULL) {
				MPI_Cancel(&req[k]);
				MPI_Wait(&req[k], MPI_STATUS_IGNORE);
			}
	}
};


/******************************* control channel ******************************/

/*
 * Layout of a solution (TAG_SOLUTION): the mixing-function index and the two
 * colliding points.  RoundState::golden has exactly this layout and is sent as is.
 * SOL_NWORDS closes the enum, so the message size is the field list itself rather
 * than a constant that has to be kept in step with it.
 *
 * The other two payloads of the engine -- a progress report (TAG_REPORT) and the
 * end-of-round statistics (MPI_Reduce) -- are both N_COUNTERS words in `enum counter`
 * order (counters.hpp): the Counters array itself, as deltas since the previous
 * report for the former and as the round's total for the latter.
 */
enum solution_field {
	SOL_I = 0, SOL_X0, SOL_X1,
	SOL_NWORDS
};

/*
 * The node's end of the control channel, on every rank.  Three tags, one per message
 * kind, so a message is told apart by its envelope and never by its length:
 *
 *     TAG_END_ROUND    rank 0 --> node    zero-length: "the round is over"
 *     TAG_REPORT       node --> rank 0    N_COUNTERS words (enum counter)
 *     TAG_SOLUTION     node --> rank 0    SOL_NWORDS words (solution_field)
 *
 * This class posts the one receive a node needs, for the end-of-round signal; it only
 * ever comes from rank 0, so the receive names its source and there is no wildcard at
 * all.  The inbound side of rank 0 -- the report and solution receives -- belongs to
 * the Controller, the only thing that ever reads them.  Rank 0 talks to itself through
 * MPI like any other node.
 *
 * Sends are plain MPI_Bsend at the call sites: the messages are short, buffered mode
 * completes locally, so there is nothing to track and the comm thread never blocks.
 * The attached buffer is per process and serves every Bsend of the engine: reports,
 * solutions, the controller's end-of-round signals and the end-of-round sentinels.
 */
class ControlChannel {
	MPI_Comm comm;
	MPI_Request req = MPI_REQUEST_NULL;
	std::vector<char> bsend_buf;

public:
	ControlChannel(const Parameters &params) : comm(params.mpi_comm)
	{
		/* Bounded in flight at once: an end-of-round signal to every node (rank 0) and
		   one sentinel per node -- OutBuffers::send_sentinels draws on this same
		   per-process buffer.  Both are zero-length.  Our own reports and solution have
		   no hard bound (reports are one-way, nothing throttles them but their pacing
		   rule); in practice a report's slot is reclaimed by the sender's own progress
		   long before the next one, so bsend_slack covers them.  A report is the largest
		   message, so sizing every slot for a full report is already generous.  A full
		   buffer is a fatal MPI_ERR_BUFFER, never a hang. */
		static_assert((int) N_COUNTERS >= (int) SOL_NWORDS, "the Bsend slots are sized for a report");
		size_t slots = 2 * (size_t) params.n_nodes + params.bsend_slack;
		size_t msg = N_COUNTERS * sizeof(u64) + MPI_BSEND_OVERHEAD;
		bsend_buf.resize(slots * msg);
		MPI_Buffer_attach(bsend_buf.data(), (int) bsend_buf.size());

		MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, comm, &req);
	}

	/* true == the end-of-round signal arrived and the receive has been reposted */
	bool poll()
	{
		int flag = 0;
		MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
		if (!flag)
			return false;
		MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, comm, &req);
		return true;
	}

	void shutdown()
	{
		if (req != MPI_REQUEST_NULL) {
			MPI_Cancel(&req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}
		void *ptr = NULL;
		int sz = 0;
		MPI_Buffer_detach(&ptr, &sz);
	}
};

}
#endif
