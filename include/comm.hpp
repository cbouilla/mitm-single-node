#ifndef MITM_COMM
#define MITM_COMM

#include <mpi.h>
#include <atomic>
#include <mutex>
#include <vector>
#include <array>
#include <memory>

#include "counters.hpp"
#include "dict.hpp"
#include "parameters.hpp"
#include "spsc.hpp"

namespace mitm {

/* spin-wait hint: these threads are pinned and own their core, so we spin rather
   than yield, but politely. */
static inline void cpu_relax()
{
#if defined(__x86_64__) || defined(__i386__)
	__builtin_ia32_pause();
#elif defined(__aarch64__)
	__asm__ __volatile__("yield");
#endif
}

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
 * the walk.  The emptiness probe is a relaxed atomic so walkers stay off the lock
 * on the common path.
 */
class CollisionQueue {
	std::mutex mtx;
	std::vector<CollisionCandidate> buf;
	size_t head = 0, tail = 0;

public:
	std::atomic<size_t> count;      /* probed without the lock on the hot path */
	std::atomic<u64> n_dropped;

	CollisionQueue(size_t capacity) : buf(capacity < 1 ? 1 : capacity), count(0), n_dropped(0) {}


	bool push(const CollisionCandidate &c)
	{
		std::lock_guard<std::mutex> lock(mtx);
		if (count.load(std::memory_order_relaxed) == buf.size()) {
			n_dropped.fetch_add(1, std::memory_order_relaxed);
			return false;
		}
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
	u64 golden[3];              /* i, x0, x1 */

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
 * Per-thread state that the comm thread also touches -- and nothing else.  Whatever a
 * walker or inserter keeps to itself (its problem wrapper, its chain state, its
 * queue, its dictionary shard) is a local of the thread function
 * or an argument to it, so what is left here is exactly the shared surface.
 */
struct alignas(64) ThreadContext {
	int role = WALKER;              /* set once at startup; the comm thread dispatches on it */
	int cpu = -1;                   /* the thread pins itself, rank 0 prints the map */

	/* see thread_state above: the comm thread asks, the thread answers */
	std::atomic<int> state;

	/* the thread tallies into this; the comm thread merges and clears it at round end */
	Counters ctr;

	/* published by the thread, read by the comm thread for its periodic report */
	std::atomic<u64> n_dp;             /* walker:   distinguished points found */
	std::atomic<u64> n_probe;          /* inserter: dictionary probes retired */
	std::atomic<u64> n_eval;           /* walker:   f evaluations, published at round end */
	std::atomic<u64> n_drop_walkerq;   /* walker:   its queue to the comm thread was full */
	std::atomic<u64> n_drop_coll;      /* inserter: the collision queue was full */

	ThreadContext()
		: state(RUNNING), n_dp(0), n_probe(0), n_eval(0), n_drop_walkerq(0), n_drop_coll(0)
	{}
};


/*************************** outgoing bulk DP buffers *************************/

/*
 * Double buffering, one pair per destination node: `ready` accumulates, `outgoing`
 * is in flight.  Never blocks -- if the previous send has not finished when `ready`
 * fills up, the point is dropped, which costs work but never correctness.
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
	u64 n_dropped = 0;
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
		if (ready[dst].size() + DP_WORDS > cap && not rotate(dst)) {
			n_dropped += 1;
			return false;
		}
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
 * itself (PcsNode::poll_incoming), because handing this class a scatter callback
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
 * Layout of a node report.  REP_NWORDS closes the enum, so the message size is the
 * field list itself rather than a constant that has to be kept in step with it.
 */
enum report_field {
	REP_NDP = 0, REP_NEVAL, REP_DROP_WALKERQ, REP_DROP_OUT, REP_DROP_INSERTERQ,
	REP_DROP_COLL, REP_GOLDEN, REP_I, REP_X0, REP_X1, REP_NPROBE,
	REP_NWORDS
};

/*
 * Layout of the end-of-round statistics.  Every thread's Counters are merged into one
 * per node, packed in this order and MPI_SUM-reduced onto rank 0 (the HyperLogLog
 * registers travel separately, under MPI_MAX).  Same trick as above: the enum closes
 * with its own length.
 */
enum round_stat {
	ST_NEVAL = 0, ST_NPOINTS_TRAILS, ST_NCOLL, ST_LEN_MIN, ST_LEN_MAX,
	ST_BAD_PROBE, ST_BAD_ROBINHOOD, ST_BAD_NONCOLLIDING, ST_BAD_COLLISION, ST_BAD_DP,
	ST_NWORDS
};

/*
 * ONE always-posted receive carries the whole control channel, in both directions
 * and on every rank: assignments from the controller to a node, and reports from a
 * node to the controller.  They are told apart by length -- an assignment is a
 * single u64, a report is REP_NWORDS of them -- so no second receive, no second tag,
 * and no separate path for rank 0 talking to itself.
 *
 * Sends are plain MPI_Bsend at the call sites: the messages are short, buffered mode
 * completes locally, so there is nothing to track and the comm thread never blocks.
 */
class ControlChannel {
	MPI_Comm comm;
	MPI_Request req = MPI_REQUEST_NULL;
	u64 buf[REP_NWORDS];
	std::vector<char> bsend_buf;

public:
	ControlChannel(const Parameters &params) : comm(params.world_comm)
	{
		static_assert(REP_NWORDS > 1, "a report must be distinguishable from an assignment by length");

		/* Worst case in flight at once: an assignment to every node (rank 0), our own
		   report, and one end-of-round sentinel per node -- OutBuffers::send_sentinels
		   draws on this same per-process buffer.  Sentinels are zero-length, so sizing
		   every slot for a full report is already generous. */
		size_t slots = 3 * (size_t) params.n_nodes + params.bsend_slack;
		size_t msg = REP_NWORDS * sizeof(u64) + MPI_BSEND_OVERHEAD;
		bsend_buf.resize(slots * msg);
		MPI_Buffer_attach(bsend_buf.data(), (int) bsend_buf.size());

		MPI_Irecv(buf, REP_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_CONTROL, comm, &req);
	}

	/*
	 * Returns how many u64 arrived, and reposts: 0 == nothing yet, 1 == an assignment
	 * (in out[0]), REP_NWORDS == a report from node *src.
	 */
	int poll(u64 out[REP_NWORDS], int *src)
	{
		int flag = 0;
		MPI_Status st;
		MPI_Test(&req, &flag, &st);
		if (!flag)
			return 0;

		int count = 0;
		MPI_Get_count(&st, MPI_UINT64_T, &count);
		for (int k = 0; k < count; k++)
			out[k] = buf[k];
		*src = st.MPI_SOURCE;

		MPI_Irecv(buf, REP_NWORDS, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_CONTROL, comm, &req);
		return count;
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
