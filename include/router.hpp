#ifndef MITM_ROUTER
#define MITM_ROUTER

#include <mpi.h>
#include <omp.h>
#include <atomic>
#include <vector>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <inttypes.h>
#include <cmath>
#include <err.h>

#include "tools.hpp"
#include "parameters.hpp"

namespace mitm {

/*
 * Route points from sender threads to globally numbered receiver threads across MPI ranks.  Procedural
 * API around one RAII handle, the thread's own Router_thread.
 *
 *   Router_thread rt = Router_Init(role, comm, tag, lossy, opts)
 *   int Router_local_rank(rt)                     rank among the category of threads (sender, receiver) of this node
 *   int Router_rank(rt)                           rank among the category of threads (sender, receiver) of all nodes
 *   int Router_local_num_send(rt)                 number of senders on this node
 *   int Router_local_num_recv(rt)                 number of receivers on this node
 *   int Router_num_send(rt)                       global number of senders
 *   int Router_num_recv(rt)                       global number of receivers
 *   Router_Push(a, b, dest, rt)                   sender: the point (a, b) to receiver `dest`
 *   Router_Close(rt)                              sender: announce "nothing more until Reset"
 *   Router_Pop(buf, n, rt)                        receiver: up to n points into buf, in order; returns how many
 *   Router_Test_drained(rt)                       receiver: nothing can arrive any more and nothing is left
 *   Router_Progress(rt)                           service: one bounded, non-blocking turn, in a loop
 *   Router_Test_quiescent(rt)                     service: the node has sent and received everything of the round
 *   Router_Stats(stats, rt)                       service: the node's tallies, exact once quiescent
 *   Router_Reset(rt)                              service, collective (an MPI_Barrier): a fresh round; no other
 *                                                 thread of the node inside a Router call meanwhile
 *
 * The round, per node: senders push then close; receivers pop until drained; the service calls Progress until
 * quiescent; a barrier; Stats, Reset; a barrier.  Lossy: no call ever waits, a point that cannot be routed at
 * once is dropped and counted.  Lossless: a sender may block, no point is ever lost.
 */

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
	alignas(64) std::atomic<size_t> head;      /* written by the consumer only */
	size_t cached_tail;                        /* consumer-private copy of tail */
	alignas(64) std::atomic<size_t> tail;      /* written by the producer only */
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

	/* producer side.  Returns false if the ring is full (the caller drops the item). */
	bool push(const Point &x)
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

	/* producer side.  The whole run or nothing, with ONE tail.store: false == not enough room. */
	bool push_all(const Point *in, size_t n)
	{
		size_t t = tail.load(std::memory_order_relaxed);
		size_t free = capacity - (t - cached_head);
		if (free < n) {
			cached_head = head.load(std::memory_order_acquire);
			free = capacity - (t - cached_head);
			if (free < n)
				return false;
		}
		for (size_t i = 0; i < n; i++)
			buf[(t + i) & mask] = in[i];
		tail.store(t + n, std::memory_order_release);
		return true;
	}

	/* consumer side */
	bool pop(Point &x)
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
	size_t pop_bulk(Point *out, size_t max)
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


/****************************** the interface's constants ******************************/

enum router_role { ROUTER_SERVICE = 0, ROUTER_SENDER = 1, ROUTER_RECEIVER = 2 };

/* Router_Stats layout */
enum router_stat {
	ROUTER_PUSHED, ROUTER_POPPED, ROUTER_SENT, ROUTER_RECV, ROUTER_LOCAL,
	ROUTER_DROPPED_SERVICE, ROUTER_DROPPED_NET, ROUTER_DROPPED_RECV,
	ROUTER_MSGS_SENT, ROUTER_MSGS_RECV, ROUTER_BYTES_SENT, ROUTER_BYTES_RECV, ROUTER_BLOCKS,
	ROUTER_STALL_OUT, ROUTER_STALL_IN, ROUTER_TURNS, ROUTER_IDLE_TURNS,
	ROUTER_STATS_SIZE
};

/* every field has a working default; every node must supply the same values */
struct Router_Opts {
	size_t block_points = 4096;         /* points per block, a power of two >= 4: one block is one message */
	size_t swc_linesize = 0;            /* points per private write-combining line, 4..block_points; 0 == auto */
	int n_recv = 32;                    /* always-posted MPI_Irecv(ANY_SOURCE) slots, and as many send slots */
	int inbox_blocks = 64;              /* a receiver's inbox, in blocks */
	int sweep_blocks = 256;             /* sealed blocks the service handles per Router_Progress: the turn's bound */
	int credit = 4;                     /* blocks in flight to one peer at most: a stalled peer holds no more send slots */
	int dests_per_node = 0;             /* destinations per node, the bench's virtual fan-out; 0 == the receivers */
	bool verbose = true;                /* rank 0 prints the sizes once connected */
};

/* Router_Init's `opts` for the defaults */
static constexpr const Router_Opts *ROUTER_DEFAULT_OPTS = NULL;


/****************************** internal constants and wire format ******************************/

static constexpr u32 ROUTER_NONE = 0xffffffffu;   /* "no block": a bare destination, an empty list, an idle slot */
static constexpr size_t ROUTER_HDR_BYTES = 64;    /* a block's header slot: the message header, one cache line */
static constexpr size_t ROUTER_DEST_WORDS = 8;    /* u64 words per destination: its two words fill one cache line */
static constexpr u32 ROUTER_MAX_L = 4096;         /* lines per block at most */
static constexpr u32 ROUTER_BATCH = 32;           /* blocks the service's stash trades with the free stack at a time */
static constexpr u64 ROUTER_LINK_BARE = 1ull << 63; /* in a sealed block's link: its destination has no block */

enum router_kind { ROUTER_DATA = 0, ROUTER_END = 1 };

/* every message starts with this, in the block's header slot; a DATA message continues with `count` points */
struct RouterMsgHdr {
	u32 kind;             /* ROUTER_DATA or ROUTER_END */
	u32 round;            /* the sender's round: an assertion, never a decision */
	u64 seq;              /* DATA: index of this message to this node; END: how many DATA messages went */
	u32 dest;             /* DATA: the receiver's local index on the target node */
	u32 count;            /* DATA: points that follow */
	u64 pad;              /* to 32 bytes */
};

static_assert(sizeof(RouterMsgHdr) == 32, "the wire header is 32 bytes");
static_assert(sizeof(Point) == 16, "a point is two words");

/* 64-byte aligned, zeroed: the zero-fill is the first touch */
static inline void *router_alloc(size_t bytes)
{
	size_t rounded = (bytes + 63) & ~(size_t) 63;
	if (rounded == 0)
		rounded = 64;
	void *p = aligned_alloc(64, rounded);
	if (p == NULL)
		err(1, "Router: aligned_alloc(%zu)", rounded);
	memset(p, 0, rounded);
	return p;
}


/****************************** the per-thread object ******************************/

class Router_node;

/*
 * A thread's object, the caller's own: Router_Init builds it in place on the thread's stack (its first touch)
 * and the caller keeps it for the team's whole life; what its calls read, with no lookup and no indirection.
 * Thread 0's owns the node.  A sender's fields are unused by a receiver and vice versa; the constructor knows
 * what to build.  A sender's private line keeps its own count in its last u64: with n < swc_linesize points in
 * slots 0..n-1 that word is n, and the push that completes the line overwrites it with the last point's val, n
 * being in a register by then.  The first cache line of the object is all a push that does not stage touches.
 * A receiver holds whole blocks: its inbox names them, its cursor reads the current one in place, and it pushes
 * a block read through onto the node's free stack.
 */
struct alignas(64) Router_thread {
	/* the push's fast path */
	Point *swc;                          /* a sender's F private lines of swc_linesize points, 64-byte aligned */
	const size_t swc_linesize;           /* the node's, copied: a push that does not stage never reads the node */
	u64 ctr[ROUTER_STATS_SIZE] = {};     /* a sender's PUSHED and DROPPED_SERVICE, a receiver's POPPED and BLOCKS; plain */
	/* identity */
	const int role;                      /* ROUTER_SERVICE, ROUTER_SENDER or ROUTER_RECEIVER */
	const int index;                     /* among the threads of its role */
	const int global_id;                 /* a sender: rank * S + index; a receiver: rank * R + index, its `dest` */
	Router_node &node;                   /* the node it belongs to */
	/* a sender's flag, on a line of its own: the service reads it every turn */
	alignas(64) std::atomic<u32> closed; /* release-stored by Router_Close, acquired by the service */
	/* a receiver's */
	RouterRing inbox;                    /* from the service: (block, count), blocks it now holds */
	u32 cur_blk = ROUTER_NONE;           /* the block being read, or NONE */
	u32 cur_off = 0;                     /* points of it delivered */
	u32 cur_count = 0;                   /* points in it */

	Router_thread(int role, int index, Router_node &rn);
	~Router_thread();
	Router_thread(const Router_thread &) = delete;
	Router_thread &operator=(const Router_thread &) = delete;
};

/* the service thread's output phase */
enum router_phase { ROUTER_OPEN, ROUTER_COLLECTING, ROUTER_FSCAN, ROUTER_FLUSHING, ROUTER_OUT_CLOSED };


/****************************** the object ******************************/

class Router_node {
public:
	/* fixed at construction */
	MPI_Comm comm;                       /* the caller's, used as is */
	const int tag;                       /* the caller's: every message of ours carries it, nothing else may */
	int rank;                            /* on it */
	int n_nodes;                         /* its size */
	const bool lossy;                    /* drop, or wait */
	Router_Opts opt;                     /* as agreed with the other nodes */

	/* fixed at connection */
	bool connected = false;              /* Router_Init has run */
	int S = 0;                           /* senders per node */
	int R = 0;                           /* receivers per node */
	int per_node = 0;                    /* destinations per node: R, or opt.dests_per_node */
	int F = 0;                           /* destinations: per_node * n_nodes */
	size_t swc_linesize = 0;             /* points per private line, resolved from the fan-out */
	u32 L = 0;                           /* lines per block */
	u32 nv_stride = 0;                   /* bytes of n_valid per block: L rounded up to a cache line */
	size_t block_bytes = 0;              /* a block: the header slot then block_points points; the longest message */
	u32 n_blocks = 0;                    /* the pool */
	size_t partial_cap = 0;              /* points per closing buffer: S * (swc_linesize - 1), whole cache lines */
	std::vector<int> roles;              /* per OpenMP thread id, its role: what connect() counts */
	std::vector<Router_thread *> senders;                     /* the S senders' objects (on their stacks), by index */
	std::vector<Router_thread *> receivers;                   /* the R receivers' objects (on their stacks), by index */
	/* a destination's line, alone on a cache line.  Word 0: (lines reserved << 32) | block id, the counter in the
	 * high word so that a runaway counter carries out of the word, not into the block id.  Word 1: the points in
	 * its closing buffer, a closer's fetch_add taking its offset. */
	std::atomic<u64> *dest = NULL;       /* F * ROUTER_DEST_WORDS */
	char *pool = NULL;                   /* n_blocks blocks of block_bytes: staging, send and receive memory alike */
	std::atomic<u8> *n_valid = NULL;     /* n_blocks * nv_stride: byte k is 0 until line k is written, then 1 */
	Point *partial = NULL;               /* F closing buffers of partial_cap points each, filled at Router_Close */
	std::atomic<u64> *blk_link = NULL;   /* per block: the next of its list, low word; sealed: its destination above */

	/* the free stack, its word alone on a cache line: (generation << 32) | head, a CAS on a stale view fails */
	alignas(64) std::atomic<u64> free_top{(u64) ROUTER_NONE};
	u64 free_pad[7] = {};                /* the rest of that line */

	/* the sealed stack, its word alone on a cache line: the sealers push, the service takes it whole */
	alignas(64) std::atomic<u32> sealed_top{ROUTER_NONE};
	u32 sealed_pad[15] = {};             /* the rest of that line */

	/* the service thread's own */
	u64 ctr[ROUTER_STATS_SIZE] = {};     /* its share of the tallies: the network, the drops it makes, the turns */
	std::vector<u32> free_list;          /* its stash: its own frees, for its own needs; a batch to or from the stack */
	u32 todo = ROUTER_NONE;              /* the service's take off the sealed stack, what is left to handle of it */
	std::vector<u64> pending;            /* popped blocks with a line still being written: dest << 32 | blk */
	std::vector<int> repairs;            /* destinations naming ROUTER_NONE, waiting for a free block */
	std::vector<int> blk_dest;           /* per block: its destination, while parked */
	std::vector<u32> blk_count;          /* per block: its points, while parked */
	std::vector<u32> park_head;          /* per receiver, then per peer: the parked FIFOs */
	std::vector<u32> park_tail;          /* their last blocks */
	std::vector<u64> inbox_pushed;       /* per receiver: blocks handed to it this round */

	std::vector<MPI_Request> out_req;    /* n_recv send slots, MPI_REQUEST_NULL when free */
	std::vector<u32> out_blk;            /* per slot: the block in flight */
	std::vector<int> out_peer;           /* per slot: its peer */
	std::vector<int> out_free;           /* the free slots */
	std::vector<int> credit_used;        /* per peer: blocks in flight to it, at most opt.credit */
	std::vector<u64> seq_sent;           /* DATA messages sent to each peer this round */
	std::vector<RouterMsgHdr> end_msg;   /* the END message to each peer, alive while its send is in flight */
	std::vector<MPI_Request> end_req;    /* its request */
	std::vector<char> end_sent;          /* issued this round */

	std::vector<MPI_Request> in_req;     /* n_recv receive slots, MPI_REQUEST_NULL when idle */
	std::vector<u32> in_blk;             /* per slot: the block posted, NONE while idle for want of a free block */
	std::vector<int> mpi_idx;            /* MPI_Testsome output */
	std::vector<MPI_Status> mpi_st;      /* MPI_Testsome output */
	std::vector<u64> n_data_recv;        /* DATA messages received, per peer */
	std::vector<long long> end_seq;      /* the END's seq, -1 until it arrives */

	/* the round */
	u32 round = 0;                       /* from 0, +1 at every Router_Reset */
	int phase = ROUTER_OPEN;             /* the output side's closure phase */
	int flush_d = 0;                     /* the F-scan's cursor: the destination */
	bool flush_built = false;            /* the fields below describe it */
	u32 flush_blk = ROUTER_NONE;         /* its installed block */
	u32 flush_k = 0;                     /* full lines of that block still to ship */
	bool flush_install = false;          /* that block went: a fresh one has to be installed */
	u32 flush_len = 0;                   /* points in its closing buffer */
	u32 flush_off = 0;                   /* of which delivered */
	std::atomic<u32> input_closed;       /* nothing can enter an inbox any more; receivers acquire it */
	bool quiescent = false;              /* what Router_Test_quiescent returns */

	Router_node(int n_threads, MPI_Comm mpi_comm, int tag, bool lossy, const Router_Opts *options);
	~Router_node();
	Router_node(const Router_node &) = delete;
	Router_node &operator=(const Router_node &) = delete;

	/* connection: Router_Init, thread 0 */
	void connect();
	void banner() const;

	/* the free stack: any thread.  Router_Push pops, Router_Pop pushes, the service trades batches with its stash */
	u32 free_pop();
	void free_push(u32 blk);
	void free_push_chain(u32 first, u32 last);
	u32 free_pop_many(u32 k, u32 *out);

	/* the sender path: Router_Push, by the push that fills a private line */
	void stage_line(Router_thread &s, int d, const Point *line);
	u32 take_fresh();
	void zero_valid(u32 blk);
	void seal(int d, u32 blk, bool bare);

	/* the service loop: Router_Progress; Router_Init and Router_Reset for the stash, repairs and receive slots */
	u32 stash_pop();
	void spill(size_t n);
	bool complete(u32 blk);
	bool place(int d, u32 blk, u32 count);
	void dispatch(int d, u32 blk, u32 count);
	void release(u32 blk);
	void park(int d, u32 blk, u32 count);
	void retry_parked(int target);
	void handle_block(int d, u32 blk);
	void repair(int d);
	void poll_out();
	void repost(int k);
	void repost_idle();
	void poll_in();
	void sweep();
	void closure();
	void start_flush(int d);
	bool fscan();
	bool flush_all();
	void check_input_closed();
	void check_quiescent();
};

static inline bool router_pow2(size_t x)
{
	return x != 0 && (x & (x - 1)) == 0;
}

/*
 * A sender gets its lines, a receiver its inbox, the service thread neither (an unused inbox has one slot).
 * The object registers itself with the node, then holds the team's last barrier: this is Router_Init's last
 * stage, and the constructor runs nowhere else.  Past the barrier the service may read every object.
 */
inline Router_thread::Router_thread(int role_, int index_, Router_node &rn)
	: swc(role_ == ROUTER_SENDER ? (Point *) router_alloc((size_t) rn.F * rn.swc_linesize * sizeof(Point)) : NULL),
	  swc_linesize(rn.swc_linesize), role(role_), index(index_),
	  global_id(role_ == ROUTER_SENDER ? rn.rank * rn.S + index_
	            : role_ == ROUTER_RECEIVER ? rn.rank * rn.R + index_ : -1),
	  node(rn), closed(0),
	  inbox(role_ == ROUTER_RECEIVER ? (size_t) rn.opt.inbox_blocks : 0)
{
	if (role == ROUTER_SENDER)
		rn.senders[index] = this;
	else if (role == ROUTER_RECEIVER)
		rn.receivers[index] = this;
	#pragma omp barrier
}

/* Thread 0's object owns the node.  The team's objects die in any order at the end of the region: a worker's
 * never touches the node here, and the node never touches a worker's. */
inline Router_thread::~Router_thread()
{
	::free(swc);
	if (role == ROUTER_SERVICE)
		delete &node;
}

/*
 * The shell, no MPI: the caller's communicator as is, `tag` for every message of ours -- the caller keeps both
 * and sends nothing else with that tag --, `lossy`, the options (NULL: the defaults) and room for the team's
 * roles.  Reached from Router_Init's single, on any thread of the team; connect() does the rest on thread 0.
 */
inline Router_node::Router_node(int n_threads, MPI_Comm mpi_comm, int tag_, bool lossy_, const Router_Opts *options)
	: comm(mpi_comm), tag(tag_), rank(-1), n_nodes(0), lossy(lossy_), opt(options ? *options : Router_Opts()),
	  input_closed(0)
{
	roles.assign(n_threads, -1);
}

/* Reached from thread 0's ~Router_thread, at the end of the region: only when nothing of ours is in flight, i.e.
 * quiescent or never used.  The workers' objects may be gone by now: nothing of theirs is read. */
inline Router_node::~Router_node()
{
	bool idle = phase == ROUTER_OPEN && ctr[ROUTER_BLOCKS] == 0 && ctr[ROUTER_MSGS_SENT] == 0;
	if (connected && not quiescent && not idle)
		errx(1, "Router: destroyed while not quiescent");
	for (size_t k = 0; k < in_req.size(); k++)
		if (in_req[k] != MPI_REQUEST_NULL) {
			MPI_Cancel(&in_req[k]);
			MPI_Wait(&in_req[k], MPI_STATUS_IGNORE);
		}
	free(dest);
	free(pool);
	free(n_valid);
	free(partial);
	free(blk_link);
}


/******************************** connection ********************************/

/*
 * Reached from Router_Init, thread 0 only, once every thread has declared its role and before any object is
 * built: the MPI checks, the agreement with the peers, the pool and the shared area, one block per destination,
 * the receives posted, a stash kept, every other block onto the free stack as one chain.
 */
inline void Router_node::connect()
{
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "Router: MPI_THREAD_FUNNELED is required (got %d)", provided);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &n_nodes);
	S = 0;
	R = 0;
	for (size_t t = 1; t < roles.size(); t++) {
		if (roles[t] == ROUTER_SENDER)
			S += 1;
		else
			R += 1;
	}
	if (S < 1 || R < 1)
		errx(1, "Router: a node needs at least one sender and one receiver (S=%d, R=%d)", S, R);

	/* every node must agree: lossy, S or R differing is fatal, an option differing means the defaults */
	const int N_AGREE = 10;
	double mine[N_AGREE] = {(double) lossy, (double) S, (double) R, (double) opt.block_points,
							(double) opt.swc_linesize, (double) opt.n_recv, (double) opt.inbox_blocks,
							(double) opt.sweep_blocks, (double) opt.credit, (double) opt.dests_per_node};
	double lo[N_AGREE], hi[N_AGREE];
	MPI_Allreduce(mine, lo, N_AGREE, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(mine, hi, N_AGREE, MPI_DOUBLE, MPI_MAX, comm);
	if (lo[0] != hi[0] || lo[1] != hi[1] || lo[2] != hi[2]) {
		if (rank == 0)
			warnx("Router: the nodes disagree on `lossy`, or have different numbers of senders or receivers");
		MPI_Abort(comm, 1);
	}
	for (int k = 3; k < N_AGREE; k++)
		if (lo[k] != hi[k]) {
			if (rank == 0)
				warnx("Router: the nodes disagree on the options; using the defaults everywhere");
			bool verbose = opt.verbose;
			opt = Router_Opts();
			opt.verbose = verbose;
			break;
		}

	if (not router_pow2(opt.block_points) || opt.block_points < 4)
		errx(1, "Router: block_points must be a power of two >= 4");
	if (opt.swc_linesize != 0 && (not router_pow2(opt.swc_linesize) || opt.swc_linesize < 4
								  || opt.swc_linesize > opt.block_points))
		errx(1, "Router: swc_linesize must be 0 or a power of two, 4 <= swc_linesize <= block_points");
	if (opt.swc_linesize != 0 && opt.block_points / opt.swc_linesize > ROUTER_MAX_L)
		errx(1, "Router: block_points / swc_linesize must be <= %u", ROUTER_MAX_L);
	if (opt.n_recv < 1 || opt.inbox_blocks < 1 || opt.sweep_blocks < 1 || opt.credit < 1)
		errx(1, "Router: n_recv, inbox_blocks, sweep_blocks and credit must be >= 1");
	if (opt.dests_per_node < 0)
		errx(1, "Router: dests_per_node must be >= 0");

	int P = n_nodes;
	credit_used.assign(P, 0);
	seq_sent.assign(P, 0);
	end_msg.assign(P, RouterMsgHdr());
	end_req.assign(P, MPI_REQUEST_NULL);
	end_sent.assign(P, 0);
	n_data_recv.assign(P, 0);
	end_seq.assign(P, -1);

	per_node = opt.dests_per_node ? opt.dests_per_node : R;
	if (per_node < R)
		errx(1, "Router: dests_per_node (%d) must be at least the receivers per node (%d)", per_node, R);
	F = per_node * n_nodes;

	size_t swc = opt.swc_linesize;
	if (swc == 0) {                       /* the lines fit L2, and a block switch is at most one line in 8 */
		size_t budget = (512 << 10) / ((size_t) F * sizeof(Point));
		swc = 4;
		while (swc * 2 <= budget && swc * 2 <= opt.block_points / 8)
			swc *= 2;
	}
	while (opt.block_points / swc > ROUTER_MAX_L)
		swc *= 2;
	swc_linesize = swc;
	L = (u32) (opt.block_points / swc);
	nv_stride = (L + 63) & ~(u32) 63;
	block_bytes = ROUTER_HDR_BYTES + opt.block_points * sizeof(Point);

	/* installed, the stash, the slots, what the receivers may hold, pending, and sealed or parked */
	size_t slack = (size_t) F > 32 * (size_t) S ? (size_t) F : 32 * (size_t) S;
	size_t nb = (size_t) F + 2 * ROUTER_BATCH + 2 * (size_t) opt.n_recv
				+ (size_t) R * ((size_t) opt.inbox_blocks + 1) + (size_t) S + slack;
	if (nb >= ROUTER_NONE)
		errx(1, "Router: too many blocks (%zu)", nb);
	n_blocks = (u32) nb;
	dest = (std::atomic<u64> *) router_alloc((size_t) F * ROUTER_DEST_WORDS * sizeof(u64));
	pool = (char *) router_alloc((size_t) n_blocks * block_bytes);
	n_valid = (std::atomic<u8> *) router_alloc((size_t) n_blocks * nv_stride);
	partial_cap = ((size_t) S * (swc_linesize - 1) + 3) & ~(size_t) 3;
	partial = (Point *) router_alloc((size_t) F * partial_cap * sizeof(Point));
	blk_link = (std::atomic<u64> *) router_alloc((size_t) n_blocks * sizeof(u64));
	free_list.reserve(n_blocks);
	for (u32 b = 0; b < n_blocks; b++) {
		blk_link[b].store(ROUTER_NONE, std::memory_order_relaxed);
		free_list.push_back(b);
	}
	blk_dest.assign(n_blocks, -1);
	blk_count.assign(n_blocks, 0);
	park_head.assign(R + n_nodes, ROUTER_NONE);
	park_tail.assign(R + n_nodes, ROUTER_NONE);
	pending.reserve(S + 1);
	inbox_pushed.assign(R, 0);

	out_req.assign(opt.n_recv, MPI_REQUEST_NULL);
	out_blk.assign(opt.n_recv, ROUTER_NONE);
	out_peer.assign(opt.n_recv, -1);
	out_free.clear();
	for (int k = 0; k < opt.n_recv; k++)
		out_free.push_back(k);
	in_req.assign(opt.n_recv, MPI_REQUEST_NULL);
	in_blk.assign(opt.n_recv, ROUTER_NONE);
	mpi_idx.assign(opt.n_recv, 0);
	mpi_st.assign(opt.n_recv, MPI_Status());
	senders.assign(S, NULL);
	receivers.assign(R, NULL);

	for (int d = 0; d < F; d++)
		dest[ROUTER_DEST_WORDS * d].store((u64) stash_pop(), std::memory_order_relaxed);
	for (int k = 0; k < opt.n_recv; k++)
		repost(k);
	if (free_list.size() > ROUTER_BATCH)
		spill(free_list.size() - ROUTER_BATCH);
	connected = true;
	if (rank == 0 && opt.verbose)
		banner();
}

/* the sizes, once.  Reached from connect(), on rank 0 only, when opt.verbose. */
inline void Router_node::banner() const
{
	double blocks = (double) n_blocks * block_bytes;
	double lines = (double) S * F * swc_linesize * sizeof(Point);
	double inboxes = (double) R * opt.inbox_blocks * sizeof(Point);
	double closing = (double) F * partial_cap * sizeof(Point);
	printf("Router: %s, %d node%s, %d sender%s and %d receiver%s per node, %d destinations (%d per node)\n",
	       lossy ? "lossy" : "lossless", n_nodes, n_nodes > 1 ? "s" : "", S, S > 1 ? "s" : "",
	       R, R > 1 ? "s" : "", F, per_node);
	printf("Router: blocks of %zu points, lines of %zu (%u lines per block), %u blocks (%.1f MB), messages of %zu"
	       " bytes\n", opt.block_points, swc_linesize, L, n_blocks, blocks / 1e6, block_bytes);
	printf("Router: private lines %.1f MB (%.1f MB per sender), closing buffers %.1f MB (%zu points per destination)\n",
	       lines / 1e6, lines / 1e6 / S, closing / 1e6, partial_cap);
	printf("Router: inboxes of %d blocks (%.1f MB), %d receive and %d send slots, credit %d per peer, stash batches"
	       " of %u, %d sealed blocks per turn\n", opt.inbox_blocks, inboxes / 1e6, opt.n_recv, opt.n_recv, opt.credit,
	       ROUTER_BATCH, opt.sweep_blocks);
}

/*
 * Collective over every thread of the team and every node, the same arguments on every thread.  Thread 0 is
 * the service whatever `role` says.  Returns the thread its object, which every later call of this thread
 * takes; the caller keeps it, named, for the team's whole life.
 */
[[nodiscard]] inline Router_thread Router_Init(int role, MPI_Comm comm, int tag, bool lossy, const Router_Opts *opts)
{
	int tid = omp_get_thread_num();
	Router_node *rn = NULL;
	#pragma omp single copyprivate(rn)
	rn = new Router_node(omp_get_num_threads(), comm, tag, lossy, opts);
	if (tid == 0)
		rn->roles[0] = ROUTER_SERVICE;
	else if (role == ROUTER_SENDER || role == ROUTER_RECEIVER)
		rn->roles[tid] = role;
	else
		errx(1, "Router_Init: role must be ROUTER_SENDER or ROUTER_RECEIVER");
	#pragma omp barrier
	if (tid == 0)
		rn->connect();
	#pragma omp barrier
	int mine = rn->roles[tid];
	int index = 0;                        /* among the threads of my role: those before me */
	for (int t = 0; t < tid; t++)
		if (rn->roles[t] == mine)
			index += 1;
	return Router_thread(mine, index, *rn);   /* registers itself, then holds the team's last barrier */
}

/* Any thread.  What the team looks like, once connected: the senders and receivers of this node and of all. */
inline int Router_local_num_send(const Router_thread &rt)
{
	return rt.node.S;
}

inline int Router_local_num_recv(const Router_thread &rt)
{
	return rt.node.R;
}

inline int Router_num_send(const Router_thread &rt)
{
	return rt.node.S * rt.node.n_nodes;
}

inline int Router_num_recv(const Router_thread &rt)
{
	return rt.node.R * rt.node.n_nodes;
}

/* Sender or receiver.  The caller's rank among the threads of its role on this node, from 0. */
inline int Router_local_rank(const Router_thread &rt)
{
	if (rt.role != ROUTER_SENDER && rt.role != ROUTER_RECEIVER)
		errx(1, "Router_local_rank: not a sender or receiver thread");
	return rt.index;
}

/* Sender or receiver.  The caller's rank among the threads of its role on all nodes, from 0: the nodes' in
 * order, each node's by local rank.  A receiver's is the `dest` that reaches it. */
inline int Router_rank(const Router_thread &rt)
{
	if (rt.role != ROUTER_SENDER && rt.role != ROUTER_RECEIVER)
		errx(1, "Router_rank: not a sender or receiver thread");
	return rt.global_id;
}


/******************************** the free stack ********************************/

/*
 * A block off the free stack, or NONE.  The generation in the high word changes at every push and pop, so
 * a CAS on a stale view fails: a `blk_link` read on the way is validated by the CAS, never trusted.  The
 * CAS acquires what the pusher did to the block before pushing it.
 * Reached from Router_Push, by the sealer: the push whose line fills a block.
 */
inline u32 Router_node::free_pop()
{
	u64 w = free_top.load(std::memory_order_acquire);
	for (;;) {
		u32 head = (u32) w;
		if (head == ROUTER_NONE)
			return ROUTER_NONE;
		u32 next = (u32) blk_link[head].load(std::memory_order_relaxed);
		u64 w2 = (((w >> 32) + 1) << 32) | next;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_acquire))
			return head;
	}
}

/* a block onto the free stack; the CAS releases everything the pusher did to it.  Reached from Router_Pop, once
 * the receiver has read a block through. */
inline void Router_node::free_push(u32 blk)
{
	u64 w = free_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[blk].store((u32) w, std::memory_order_relaxed);
		u64 w2 = (((w >> 32) + 1) << 32) | blk;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_relaxed))
			return;
	}
}

/* a chain the caller linked through blk_link, `first` down to `last`, onto the free stack in one CAS.  Reached
 * from Router_Init (thread 0) and Router_Progress, when the service's stash spills. */
inline void Router_node::free_push_chain(u32 first, u32 last)
{
	u64 w = free_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[last].store((u32) w, std::memory_order_relaxed);
		u64 w2 = (((w >> 32) + 1) << 32) | first;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_relaxed))
			return;
	}
}

/*
 * Up to k blocks off the free stack in one CAS, after walking k links.  A walk over a stale view is
 * harmless: it is bounded, every link it reads is a block id or NONE, and its CAS fails.
 * Reached from Router_Init, Router_Progress and Router_Reset, whenever the service's stash runs empty.
 */
inline u32 Router_node::free_pop_many(u32 k, u32 *out)
{
	u64 w = free_top.load(std::memory_order_acquire);
	for (;;) {
		u32 n = 0;
		u32 cur = (u32) w;
		while (n < k && cur != ROUTER_NONE) {
			out[n] = cur;
			n += 1;
			cur = (u32) blk_link[cur].load(std::memory_order_relaxed);
		}
		if (n == 0)
			return 0;
		u64 w2 = (((w >> 32) + 1) << 32) | cur;
		if (free_top.compare_exchange_weak(w, w2, std::memory_order_acq_rel, std::memory_order_acquire))
			return n;
	}
}


/******************************** the sender path ********************************/

/* a free block for the sealer, off the free stack: lossless waits for one, lossy gets NONE.  Reached from
 * Router_Push, by the push whose line fills a block. */
inline u32 Router_node::take_fresh()
{
	u32 blk = free_pop();
	if (blk != ROUTER_NONE || lossy)
		return blk;
	while ((blk = free_pop()) == ROUTER_NONE)
		cpu_relax();
	return blk;
}

/* a block about to be installed: none of its lines is written.  Reached from Router_Push (the sealer's fresh
 * block), Router_Progress (a repair, the F-scan's fresh block) and Router_Reset (a repair). */
inline void Router_node::zero_valid(u32 blk)
{
	for (u32 k = 0; k < L; k++)
		n_valid[(size_t) nv_stride * blk + k].store(0, std::memory_order_relaxed);
}

/*
 * A sealed block onto the sealed stack, its destination in the link's high word with the BARE bit when the
 * sealer left the slot without a block; the CAS releases that slot's store.  No generation word: a push on a
 * stale view fails or links to the true top.  Reached from Router_Push, by the push whose line fills a block.
 */
inline void Router_node::seal(int d, u32 blk, bool bare)
{
	u64 hi = ((u64) d << 32) | (bare ? ROUTER_LINK_BARE : 0);
	u32 top = sealed_top.load(std::memory_order_relaxed);
	for (;;) {
		blk_link[blk].store(hi | top, std::memory_order_relaxed);
		if (sealed_top.compare_exchange_weak(top, blk, std::memory_order_release, std::memory_order_relaxed))
			return;
	}
}

/*
 * Reserve one line slot in the block of destination `d` and write the full line there, so that a block in
 * flight never has a hole; seal the block when the reservation lands on its end; drop the line in lossy
 * mode when the destination has no block or a sealer is installing one: lossy never waits.
 * Reached from Router_Push, by the push that fills the sender's private line for `d`.
 */
inline void Router_node::stage_line(Router_thread &s, int d, const Point *line)
{
	std::atomic<u64> &next = dest[ROUTER_DEST_WORDS * d];
	int n = (int) swc_linesize;
	for (;;) {
		u64 v = next.fetch_add(1ull << 32, std::memory_order_acq_rel);
		u32 k = (u32) (v >> 32);
		u32 blk = (u32) v;
		if (blk == ROUTER_NONE) {         /* the service thread repairs it, once it pops the note */
			s.ctr[ROUTER_DROPPED_SERVICE] += n;
			return;
		}
		if (k < L) {
			Point *pts = (Point *) (pool + (size_t) blk * block_bytes + ROUTER_HDR_BYTES);
			memcpy(pts + (size_t) k * swc_linesize, line, (size_t) n * sizeof(Point));
			n_valid[(size_t) nv_stride * blk + k].store(1, std::memory_order_release);
			return;
		}
		if (k == L) {
			u32 fresh = take_fresh();
			if (fresh != ROUTER_NONE)
				zero_valid(fresh);
			next.store((u64) fresh, std::memory_order_release);
			seal(d, blk, fresh == ROUTER_NONE);
			continue;
		}
		if (lossy) {                      /* k > L: a sealer is installing, and lossy waits for nobody */
			s.ctr[ROUTER_DROPPED_SERVICE] += n;
			return;
		}
		for (;;) {                        /* at exactly L nobody is installing, and I may become the sealer */
			u64 w = next.load(std::memory_order_acquire);
			if ((u32) (w >> 32) <= L)
				break;
			cpu_relax();
		}
	}
}

/* Sender.  Lossy: wait-free, may drop the line; lossless: may spin, never loses the point. */
inline void Router_Push(u64 a, u64 b, int d, Router_thread &rt)
{
	Point *line = rt.swc + (size_t) d * rt.swc_linesize;
	u64 &count = line[rt.swc_linesize - 1].val;
	u64 n = count;
	line[n].key = a;
	line[n].val = b;                  /* the line's last point lands on the count, read already */
	n += 1;
	rt.ctr[ROUTER_PUSHED] += 1;
	if (n < rt.swc_linesize) {
		count = n;
		return;
	}
	rt.node.stage_line(rt, d, line);
	count = 0;                        /* after the copy: the word held a point until now */
}

/*
 * Sender.  Append every partial line's points to its destination's closing buffer, contiguous with the other
 * senders' thanks to the offset a fetch_add hands out, then announce; nothing more from this sender until
 * Router_Reset.  Nothing partial ever enters a block before the F-scan copies the closing buffers into blocks.
 */
inline void Router_Close(Router_thread &rt)
{
	Router_node &rn = rt.node;
	for (int d = 0; d < rn.F; d++) {
		Point *line = rt.swc + (size_t) d * rt.swc_linesize;
		u64 &count = line[rt.swc_linesize - 1].val;
		if (count == 0)
			continue;
		u64 off = rn.dest[ROUTER_DEST_WORDS * d + 1].fetch_add(count, std::memory_order_relaxed);
		memcpy(rn.partial + (size_t) d * rn.partial_cap + off, line, (size_t) count * sizeof(Point));
		count = 0;
	}
	rt.closed.store(1, std::memory_order_release);
}


/******************************** the receiver's calls ********************************/

/*
 * Receiver.  Up to `buffer_size` points into `buffer` (2 u64 each), out of the blocks the inbox names, in
 * order; returns how many.  A block read through goes onto the free stack, whose CAS is the release that
 * orders these reads before the block's reuse, and is counted after that CAS, which keeps the count behind it.
 */
inline size_t Router_Pop(u64 *buffer, size_t buffer_size, Router_thread &rt)
{
	Router_node &rn = rt.node;
	size_t got = 0;
	while (got < buffer_size) {
		if (rt.cur_blk == ROUTER_NONE) {
			Point e;
			if (not rt.inbox.pop(e))
				break;
			rt.cur_blk = (u32) e.key;
			rt.cur_count = (u32) e.val;
			rt.cur_off = 0;
		}
		size_t take = rt.cur_count - rt.cur_off;
		if (take > buffer_size - got)
			take = buffer_size - got;
		const Point *pts = (const Point *) (rn.pool + (size_t) rt.cur_blk * rn.block_bytes + ROUTER_HDR_BYTES);
		memcpy(buffer + 2 * got, pts + rt.cur_off, take * sizeof(Point));
		got += take;
		rt.cur_off += (u32) take;
		if (rt.cur_off == rt.cur_count) {
			rn.free_push(rt.cur_blk);
			rt.ctr[ROUTER_BLOCKS] += 1;
			rt.cur_blk = ROUTER_NONE;
		}
	}
	rt.ctr[ROUTER_POPPED] += got;
	return got;
}

/* Receiver.  Nothing can arrive any more and nothing is left to read: the flag first, then emptiness, so that
 * a block pushed between the two reads is still popped. */
inline bool Router_Test_drained(Router_thread &rt)
{
	if (rt.node.input_closed.load(std::memory_order_acquire) == 0)
		return false;
	return rt.cur_blk == ROUTER_NONE && rt.inbox.empty();
}


/******************************** the service thread ********************************/

/* a block from the stash, refilled from the free stack by the batch when empty; NONE if both are empty.  Reached
 * from Router_Init (the installs, the receives), Router_Progress (a repost, a repair, the F-scan) and
 * Router_Reset (a repair, a repost). */
inline u32 Router_node::stash_pop()
{
	if (free_list.empty()) {
		u32 got[ROUTER_BATCH];
		u32 n = free_pop_many(ROUTER_BATCH, got);
		for (u32 i = 0; i < n; i++)
			free_list.push_back(got[i]);
		if (n == 0)
			return ROUTER_NONE;
	}
	u32 blk = free_list.back();
	free_list.pop_back();
	return blk;
}

/* the top n blocks of the stash onto the free stack, linked here, pushed by one CAS.  Reached from Router_Init,
 * once the stash exceeds a batch, and Router_Progress, once a release takes it past two. */
inline void Router_node::spill(size_t n)
{
	size_t keep = free_list.size() - n;
	for (size_t i = keep; i + 1 < free_list.size(); i++)
		blk_link[free_list[i]].store(free_list[i + 1], std::memory_order_relaxed);
	free_push_chain(free_list[keep], free_list.back());
	free_list.resize(keep);
}

/*
 * Are all L lines of a sealed block written?  A zero n_valid byte means not yet, never a torn value, each
 * line being its own release/acquire pair.  A written line is full: a block in flight has no hole.
 * Reached from Router_Progress, for every sealed block the take or the pending list yields.
 */
inline bool Router_node::complete(u32 blk)
{
	for (u32 k = 0; k < L; k++)
		if (n_valid[(size_t) nv_stride * blk + k].load(std::memory_order_acquire) == 0)
			return false;
	return true;
}

/*
 * A block to its destination's owner: the local receiver's inbox, or an Isend from block memory to its node.
 * True == the block changed hands, and the caller must forget it.  False == no room: the inbox is full, or
 * no send slot is free, or the peer has its credit in flight.
 * Reached from Router_Progress, for every block dispatched and every parked block retried.
 */
inline bool Router_node::place(int d, u32 blk, u32 count)
{
	int node = d / per_node;
	int local = (d % per_node) % R;
	if (node == rank) {
		Point e = {blk, count};
		if (not receivers[local]->inbox.push(e))
			return false;
		inbox_pushed[local] += 1;
		ctr[ROUTER_LOCAL] += count;
		return true;
	}
	if (out_free.empty() || credit_used[node] >= opt.credit)
		return false;
	int k = out_free.back();
	out_free.pop_back();
	char *hdr = pool + (size_t) blk * block_bytes;
	RouterMsgHdr h = {ROUTER_DATA, round, seq_sent[node], (u32) local, count, 0};
	memcpy(hdr, &h, sizeof(h));
	int len = (int) (ROUTER_HDR_BYTES + (size_t) count * sizeof(Point));
	MPI_Isend(hdr, len, MPI_BYTE, node, tag, comm, &out_req[k]);
	out_blk[k] = blk;
	out_peer[k] = node;
	credit_used[node] += 1;
	seq_sent[node] += 1;
	ctr[ROUTER_SENT] += count;
	ctr[ROUTER_MSGS_SENT] += 1;
	ctr[ROUTER_BYTES_SENT] += len;
	return true;
}

/* a complete block to its destination: placed, parked (lossless) or dropped (lossy), never kept by the caller.
 * Reached from Router_Progress: a sealed block once complete, a DATA message received, the F-scan's blocks. */
inline void Router_node::dispatch(int d, u32 blk, u32 count)
{
	if (place(d, blk, count))
		return;
	if (lossy) {
		ctr[(d / per_node == rank) ? ROUTER_DROPPED_RECV : ROUTER_DROPPED_NET] += count;
		release(blk);
		return;
	}
	park(d, blk, count);
}

/* a block the service is done with: into its stash, which spills a batch to the free stack when it grows.  Reached
 * from Router_Progress: a send that completed, or (lossy) a block dropped for want of room. */
inline void Router_node::release(u32 blk)
{
	ctr[ROUTER_BLOCKS] += 1;
	free_list.push_back(blk);
	if (free_list.size() > 2 * ROUTER_BATCH)
		spill(ROUTER_BATCH);
}

/* lossless: keep a block whose target is full, in arrival order, until the target has room; parking rather
 * than refusing to pop is what keeps one slow receiver from stalling every other destination.  Reached from
 * Router_Progress, lossless only, when a block finds its inbox full or its peer out of slots or credit. */
inline void Router_node::park(int d, u32 blk, u32 count)
{
	int node = d / per_node;
	int target = (node == rank) ? (d % per_node) % R : R + node;
	blk_dest[blk] = d;
	blk_count[blk] = count;
	blk_link[blk].store(ROUTER_NONE, std::memory_order_relaxed);
	if (park_head[target] == ROUTER_NONE)
		park_head[target] = blk;
	else
		blk_link[park_tail[target]].store(blk, std::memory_order_relaxed);
	park_tail[target] = blk;
	ctr[ROUTER_STALL_OUT] += 1;
}

/* one target's parked blocks, in order, as far as they place.  Reached from Router_Progress, every turn for every
 * target; a no-op unless lossless with blocks parked.  careful: the link is read before the block is placed,
 * since once placed its receiver may relink it onto the free stack. */
inline void Router_node::retry_parked(int target)
{
	while (park_head[target] != ROUTER_NONE) {
		u32 blk = park_head[target];
		u32 after = (u32) blk_link[blk].load(std::memory_order_relaxed);
		if (not place(blk_dest[blk], blk, blk_count[blk]))
			return;
		park_head[target] = after;
		if (after == ROUTER_NONE)
			park_tail[target] = ROUTER_NONE;
	}
}

/* a sealed block: dispatch it, or pend it until its last line is written.  Reached from Router_Progress, for
 * every block of the take and every block pended the turn before. */
inline void Router_node::handle_block(int d, u32 blk)
{
	if (not complete(blk)) {
		pending.push_back(((u64) d << 32) | blk);
		return;
	}
	dispatch(d, blk, (u32) opt.block_points);
}

/* lossy: a destination left without a block gets one from the stash, if it still has none.  Reached from
 * Router_Progress (a sealed block flagged BARE, the repairs list, a bare slot at the F-scan) and from
 * Router_Reset (the repairs list, which must then empty). */
inline void Router_node::repair(int d)
{
	std::atomic<u64> &next = dest[ROUTER_DEST_WORDS * d];
	u64 v = next.load(std::memory_order_acquire);
	if ((u32) v != ROUTER_NONE)
		return;
	u32 fresh = stash_pop();
	if (fresh == ROUTER_NONE) {
		repairs.push_back(d);
		return;
	}
	zero_valid(fresh);
	while (not next.compare_exchange_weak(v, (u64) fresh, std::memory_order_acq_rel))
		if ((u32) v != ROUTER_NONE) {     /* a block appeared meanwhile: keep this one */
			free_list.push_back(fresh);
			return;
		}
}

/* the sends that completed: their blocks are the service's again, their slots and their peer's credit too.
 * Reached from Router_Progress, every turn, first. */
inline void Router_node::poll_out()
{
	int outcount = 0;
	MPI_Testsome(opt.n_recv, out_req.data(), &outcount, mpi_idx.data(), mpi_st.data());
	if (outcount == MPI_UNDEFINED)
		return;
	for (int t = 0; t < outcount; t++) {
		int k = mpi_idx[t];
		release(out_blk[k]);
		credit_used[out_peer[k]] -= 1;
		out_blk[k] = ROUTER_NONE;
		out_free.push_back(k);
	}
}

/* post a block from the stash on receive slot k; none == the slot idles, and repost_idle retries every turn.
 * Reached from Router_Init (every slot), Router_Progress (a slot whose DATA block went to its receiver, and
 * the idle slots) and Router_Reset (the idle slots). */
inline void Router_node::repost(int k)
{
	u32 blk = stash_pop();
	if (blk == ROUTER_NONE) {
		ctr[ROUTER_STALL_IN] += 1;
		return;
	}
	in_blk[k] = blk;
	MPI_Irecv(pool + (size_t) blk * block_bytes, (int) block_bytes, MPI_BYTE, MPI_ANY_SOURCE, tag, comm, &in_req[k]);
}

/* a block on every idle receive slot, as far as the stash and the free stack allow.  Reached from
 * Router_Progress, every turn, and Router_Reset, once. */
inline void Router_node::repost_idle()
{
	for (int k = 0; k < opt.n_recv; k++)
		if (in_blk[k] == ROUTER_NONE) {
			repost(k);
			if (in_blk[k] == ROUTER_NONE)
				return;
		}
}

/*
 * Take delivery of what arrived: an END is recorded and its block reposted; a DATA block goes to its receiver
 * as it is -- parked if the inbox is full (lossless), dropped (lossy) -- and its slot gets a fresh block.
 * Reached from Router_Progress, every turn.
 */
inline void Router_node::poll_in()
{
	int outcount = 0;
	MPI_Testsome(opt.n_recv, in_req.data(), &outcount, mpi_idx.data(), mpi_st.data());
	if (outcount == MPI_UNDEFINED)
		return;
	for (int t = 0; t < outcount; t++) {
		int k = mpi_idx[t];
		u32 blk = in_blk[k];
		int src = mpi_st[t].MPI_SOURCE;
		int len = 0;
		MPI_Get_count(&mpi_st[t], MPI_BYTE, &len);
		RouterMsgHdr h;
		memcpy(&h, pool + (size_t) blk * block_bytes, sizeof(h));
		if (h.round != round)
			errx(1, "Router: rank %d got a round %u message in round %u", rank, h.round, round);
		ctr[ROUTER_MSGS_RECV] += 1;
		ctr[ROUTER_BYTES_RECV] += len;
		if (h.kind == ROUTER_END) {
			end_seq[src] = (long long) h.seq;
			MPI_Irecv(pool + (size_t) blk * block_bytes, (int) block_bytes, MPI_BYTE, MPI_ANY_SOURCE, tag, comm,
			          &in_req[k]);
			continue;
		}
		n_data_recv[src] += 1;
		ctr[ROUTER_RECV] += h.count;
		in_blk[k] = ROUTER_NONE;
		dispatch(rank * per_node + (int) h.dest, blk, h.count);
		repost(k);
	}
}

/* the pending and repair lists, then the sealed blocks: the stack taken whole in one exchange once the previous
 * take is handled, and at most opt.sweep_blocks of the take handled per turn, newest first, so that the turn
 * stays bounded.  Reached from Router_Progress, every turn. */
inline void Router_node::sweep()
{
	if (not pending.empty()) {
		std::vector<u64> again;
		again.swap(pending);
		for (size_t i = 0; i < again.size(); i++)
			handle_block((int) (again[i] >> 32), (u32) again[i]);
	}
	if (not repairs.empty()) {
		std::vector<int> again;
		again.swap(repairs);
		for (size_t i = 0; i < again.size(); i++)
			repair(again[i]);
	}
	if (todo == ROUTER_NONE)
		todo = sealed_top.exchange(ROUTER_NONE, std::memory_order_acquire);
	for (int b = 0; b < opt.sweep_blocks && todo != ROUTER_NONE; b++) {
		u32 blk = todo;
		u64 link = blk_link[blk].load(std::memory_order_relaxed);   /* careful: handled, the block may be relinked */
		todo = (u32) link;
		int d = (int) ((link >> 32) & 0x7fffffffu);
		if (link & ROUTER_LINK_BARE)
			repair(d);
		handle_block(d, blk);
	}
}

/*
 * Destination `d` at closure, once every sender is closed: what its installed block holds, and how many
 * points its closing buffer holds.  A NONE slot gets a block from the stash if there is one.
 * Reached from Router_Progress in the F-scan phase, once per destination and round.
 */
inline void Router_node::start_flush(int d)
{
	u64 v = dest[ROUTER_DEST_WORDS * d].load(std::memory_order_acquire);
	u32 k = (u32) (v >> 32);
	u32 blk = (u32) v;
	if (blk == ROUTER_NONE) {
		repair(d);
		k = 0;
	} else {
		if (k > L)
			errx(1, "Router: destination %d has %u lines reserved at closure", d, k);
		for (u32 i = 0; i < k; i++)
			if (n_valid[(size_t) nv_stride * blk + i].load(std::memory_order_acquire) == 0)
				errx(1, "Router: destination %d has an unwritten line at closure", d);
	}
	flush_blk = blk;
	flush_k = k;
	flush_install = false;
	flush_len = (u32) dest[ROUTER_DEST_WORDS * d + 1].load(std::memory_order_relaxed);
	flush_off = 0;
	flush_built = true;
}

/*
 * The closing runs, one destination per step.  The installed block's full lines go as one short block, the
 * block itself, and a fresh one is installed in its place; then the closing buffer is copied into fresh
 * blocks, at most a block at a time, each dispatched -- the service's one copy, per round.  Both need free
 * blocks: lossless stops where there is none and resumes there next turn, lossy drops the run.
 * Reached from Router_Progress, every turn of the F-scan phase, until it returns true.
 */
inline bool Router_node::fscan()
{
	for (; flush_d < F; flush_d++) {
		if (not flush_built)
			start_flush(flush_d);
		if (flush_k > 0) {
			dispatch(flush_d, flush_blk, flush_k * (u32) swc_linesize);
			flush_k = 0;
			flush_install = true;
		}
		if (flush_install) {
			u32 fresh = stash_pop();
			if (fresh == ROUTER_NONE) {
				if (not lossy)
					return false;
				dest[ROUTER_DEST_WORDS * flush_d].store((u64) ROUTER_NONE, std::memory_order_release);
				repair(flush_d);
			} else {
				zero_valid(fresh);
				dest[ROUTER_DEST_WORDS * flush_d].store((u64) fresh, std::memory_order_release);
			}
			flush_install = false;
		}
		const Point *buf = partial + (size_t) flush_d * partial_cap;
		while (flush_off < flush_len) {
			u32 chunk = flush_len - flush_off;
			if (chunk > opt.block_points)
				chunk = (u32) opt.block_points;
			u32 blk = stash_pop();
			if (blk == ROUTER_NONE) {
				if (not lossy)
					return false;
				ctr[(flush_d / per_node == rank) ? ROUTER_DROPPED_RECV : ROUTER_DROPPED_NET] += chunk;
				flush_off += chunk;
				continue;
			}
			memcpy(pool + (size_t) blk * block_bytes + ROUTER_HDR_BYTES, buf + flush_off, (size_t) chunk * sizeof(Point));
			dispatch(flush_d, blk, chunk);
			flush_off += chunk;
		}
		dest[ROUTER_DEST_WORDS * flush_d + 1].store(0, std::memory_order_relaxed);
		flush_built = false;
	}
	return true;
}

/*
 * Send the ENDs, test every request; true once nothing of ours is in flight.  A peer's END waits for its
 * parked blocks: they still have to go out, and the END must follow them.
 * Reached from Router_Progress, every turn of the flushing phase, until it returns true.
 */
inline bool Router_node::flush_all()
{
	bool done = true;
	for (int p = 0; p < n_nodes; p++) {
		if (p == rank)
			continue;
		if (park_head[R + p] != ROUTER_NONE) {
			done = false;
			continue;
		}
		if (not end_sent[p]) {
			end_msg[p] = RouterMsgHdr{ROUTER_END, round, seq_sent[p], 0, 0, 0};
			MPI_Isend(&end_msg[p], sizeof(RouterMsgHdr), MPI_BYTE, p, tag, comm, &end_req[p]);
			end_sent[p] = 1;
			ctr[ROUTER_MSGS_SENT] += 1;
			ctr[ROUTER_BYTES_SENT] += sizeof(RouterMsgHdr);
		}
		if (end_req[p] != MPI_REQUEST_NULL) {
			int flag = 0;
			MPI_Test(&end_req[p], &flag, MPI_STATUS_IGNORE);
			if (not flag)
				done = false;
		}
	}
	if (out_free.size() != (size_t) opt.n_recv)
		done = false;
	return done;
}

/* nothing can enter an inbox any more: local flush done, no block parked on a receiver, every peer's END and
 * all its DATA in.  The END's seq is compared with the DATA count because MPI orders matching, not
 * completion, across posted receives.  Reached from Router_Progress, every turn; a no-op before the flushing
 * phase. */
inline void Router_node::check_input_closed()
{
	if (phase < ROUTER_FLUSHING)
		return;
	for (int r = 0; r < R; r++)
		if (park_head[r] != ROUTER_NONE)
			return;
	for (int p = 0; p < n_nodes; p++)
		if (p != rank && (end_seq[p] < 0 || (u64) end_seq[p] != n_data_recv[p]))
			return;
	input_closed.store(1, std::memory_order_release);
}

/* out closed, input closed, and every block handed to a receiver given back: nobody holds anything.  The
 * receiver's count is plain and read racily; it is written after the push it counts, so it is never early.
 * Reached from Router_Progress, every turn; a no-op until the output side is closed and the input closed. */
inline void Router_node::check_quiescent()
{
	if (phase != ROUTER_OUT_CLOSED || input_closed.load(std::memory_order_relaxed) == 0)
		return;
	for (int r = 0; r < R; r++)
		if (inbox_pushed[r] != receivers[r]->ctr[ROUTER_BLOCKS])
			return;
	quiescent = true;
}

/* the output side's phases, one step per turn at most.  Reached from Router_Progress, every turn; a phase
 * advances only once its condition holds, the first being every sender closed. */
inline void Router_node::closure()
{
	switch (phase) {
	case ROUTER_OPEN:
		for (int s = 0; s < S; s++)
			if (senders[s]->closed.load(std::memory_order_acquire) == 0)
				return;
		phase = ROUTER_COLLECTING;
		return;
	case ROUTER_COLLECTING:
		if (not pending.empty() || todo != ROUTER_NONE || sealed_top.load(std::memory_order_acquire) != ROUTER_NONE)
			return;
		phase = ROUTER_FSCAN;
		flush_d = 0;
		return;
	case ROUTER_FSCAN:
		if (fscan())
			phase = ROUTER_FLUSHING;
		return;
	case ROUTER_FLUSHING:
		if (flush_all())
			phase = ROUTER_OUT_CLOSED;
		return;
	default:
		return;
	}
}

/* Service.  One bounded, non-blocking turn: every part runs whatever the phase, so no stall couples two; the
 * sends that completed come first, so that this turn's receive slots have their blocks. */
inline void Router_Progress(Router_thread &rt)
{
	if (rt.role != ROUTER_SERVICE)
		errx(1, "Router_Progress: not the service thread");
	Router_node &rn = rt.node;
	if (rn.quiescent)
		return;
	u64 before = rn.ctr[ROUTER_BLOCKS] + rn.ctr[ROUTER_MSGS_SENT] + rn.ctr[ROUTER_MSGS_RECV] + rn.ctr[ROUTER_RECV];
	rn.poll_out();
	for (int t = 0; t < rn.R + rn.n_nodes; t++)
		rn.retry_parked(t);
	rn.poll_in();
	rn.repost_idle();
	rn.sweep();
	rn.closure();
	rn.check_input_closed();
	rn.check_quiescent();
	rn.ctr[ROUTER_TURNS] += 1;
	u64 after = rn.ctr[ROUTER_BLOCKS] + rn.ctr[ROUTER_MSGS_SENT] + rn.ctr[ROUTER_MSGS_RECV] + rn.ctr[ROUTER_RECV];
	if (after == before)
		rn.ctr[ROUTER_IDLE_TURNS] += 1;
}

/* Service.  Every piece of the node's state, for a run that does not end; the reads are racy on purpose. */
inline void Router_Dump(FILE *f, const Router_thread &rt)
{
	if (rt.role != ROUTER_SERVICE)
		errx(1, "Router_Dump: not the service thread");
	Router_node &rn = rt.node;
	u64 top = rn.free_top.load();
	fprintf(f, "rank %d round %u phase %d input_closed %u quiescent %d flush_d %d run %u/%u pending %zu"
	        " repairs %zu sealed top %u todo %u stash %zu stack head %u gen %u send slots free %zu\n", rn.rank,
	        rn.round, rn.phase, rn.input_closed.load(), (int) rn.quiescent, rn.flush_d, rn.flush_off, rn.flush_len,
	        rn.pending.size(), rn.repairs.size(), rn.sealed_top.load(), rn.todo, rn.free_list.size(), (u32) top,
	        (u32) (top >> 32), rn.out_free.size());
	for (int s = 0; s < rn.S; s++) {
		Router_thread &sd = *rn.senders[s];
		fprintf(f, "  sender %d: closed %u pushed %" PRIu64 " dropped %" PRIu64 "\n", s, sd.closed.load(),
		        sd.ctr[ROUTER_PUSHED], sd.ctr[ROUTER_DROPPED_SERVICE]);
	}
	for (int r = 0; r < rn.R; r++) {
		Router_thread &rr = *rn.receivers[r];
		fprintf(f, "  receiver %d: inbox %zu/%zu cur %u (%u/%u) pushed %" PRIu64 " given back %" PRIu64 " parked %u\n",
		        r, (size_t) rr.inbox.head.load(), (size_t) rr.inbox.tail.load(), rr.cur_blk, rr.cur_off, rr.cur_count,
		        rn.inbox_pushed[r], rr.ctr[ROUTER_BLOCKS], rn.park_head[r]);
	}
	for (int d = 0; d < rn.F; d++) {
		u64 v = rn.dest[ROUTER_DEST_WORDS * d].load();
		fprintf(f, "  dest %d: block %u lines %u", d, (u32) v, (u32) (v >> 32));
		if ((u32) v != ROUTER_NONE)
			for (u32 k = 0; k < rn.L && k < 64; k++)
				fprintf(f, " %u", rn.n_valid[(size_t) rn.nv_stride * (u32) v + k].load());
		fprintf(f, "\n");
	}
	for (int p = 0; p < rn.n_nodes; p++) {
		if (p == rn.rank)
			continue;
		fprintf(f, "  peer %d: seq_sent %" PRIu64 " n_data_recv %" PRIu64 " end_seq %lld end_sent %d end_req %d"
		        " credit_used %d parked %u\n", p, rn.seq_sent[p], rn.n_data_recv[p], rn.end_seq[p],
		        (int) rn.end_sent[p], rn.end_req[p] != MPI_REQUEST_NULL, rn.credit_used[p], rn.park_head[rn.R + p]);
	}
	for (int k = 0; k < rn.opt.n_recv; k++) {
		if (rn.out_blk[k] != ROUTER_NONE)
			fprintf(f, "  send slot %d: block %u to peer %d\n", k, rn.out_blk[k], rn.out_peer[k]);
		if (rn.in_blk[k] == ROUTER_NONE)
			fprintf(f, "  receive slot %d: idle\n", k);
	}
	fflush(f);
}

/* Service.  This node has sent and received everything of the round and its receivers hold nothing. */
inline bool Router_Test_quiescent(const Router_thread &rt)
{
	return rt.node.quiescent;
}

/* Service.  The node's tallies; exact once quiescent and the caller's barrier has passed, a snapshot before. */
inline void Router_Stats(u64 *stats, const Router_thread &rt)
{
	if (rt.role != ROUTER_SERVICE)
		errx(1, "Router_Stats: not the service thread");
	Router_node &rn = rt.node;
	for (int k = 0; k < ROUTER_STATS_SIZE; k++)
		stats[k] = rn.ctr[k];
	for (int s = 0; s < rn.S; s++)
		for (int k = 0; k < ROUTER_STATS_SIZE; k++)
			stats[k] += rn.senders[s]->ctr[k];
	for (int r = 0; r < rn.R; r++)
		for (int k = 0; k < ROUTER_STATS_SIZE; k++)
			stats[k] += rn.receivers[r]->ctr[k];
}

/*
 * Service, collective.  Back to a fresh round: the counters and flags, then an MPI_Barrier so that no
 * peer starts the next round before this node has left this one.  No other thread of
 * the node may be inside a Router call meanwhile.
 */
inline void Router_Reset(Router_thread &rt)
{
	if (rt.role != ROUTER_SERVICE)
		errx(1, "Router_Reset: not the service thread");
	Router_node &rn = rt.node;
	if (not rn.quiescent)
		errx(1, "Router_Reset: not quiescent");
	if (not rn.repairs.empty()) {         /* lossy: destinations left without a block at the F-scan */
		std::vector<int> again;
		again.swap(rn.repairs);
		for (size_t i = 0; i < again.size(); i++)
			rn.repair(again[i]);
		if (not rn.repairs.empty())
			errx(1, "Router_Reset: a destination has no block");
	}
	rn.repost_idle();
#ifdef ROUTER_PARANOID
	size_t accounted = rn.F + rn.free_list.size();
	for (int k = 0; k < rn.opt.n_recv; k++)
		if (rn.in_blk[k] != ROUTER_NONE)
			accounted += 1;
	size_t walked = 0;                    /* the free stack, nobody else touching it now */
	for (u32 b = (u32) rn.free_top.load(); b != ROUTER_NONE && walked <= rn.n_blocks; b = (u32) rn.blk_link[b].load())
		walked += 1;
	accounted += walked;
	if (accounted != rn.n_blocks)
		errx(1, "Router_Reset: %zu of %u blocks accounted for (%zu on the stack)", accounted, rn.n_blocks, walked);
	if (rn.out_free.size() != (size_t) rn.opt.n_recv)
		errx(1, "Router_Reset: a send is still in flight");
	for (int p = 0; p < rn.n_nodes; p++)
		if (rn.credit_used[p] != 0)
			errx(1, "Router_Reset: credit in use");
	if (not rn.pending.empty())
		errx(1, "Router_Reset: pending work left");
	if (rn.todo != ROUTER_NONE || rn.sealed_top.load() != ROUTER_NONE)
		errx(1, "Router_Reset: sealed blocks left");
	for (int d = 0; d < rn.F; d++)
		if (rn.dest[ROUTER_DEST_WORDS * d + 1].load(std::memory_order_relaxed) != 0)
			errx(1, "Router_Reset: a closing buffer was left behind");
#endif
	for (int k = 0; k < ROUTER_STATS_SIZE; k++)
		rn.ctr[k] = 0;
	for (int s = 0; s < rn.S; s++) {
		Router_thread &sd = *rn.senders[s];
		for (int k = 0; k < ROUTER_STATS_SIZE; k++)
			sd.ctr[k] = 0;
		sd.closed.store(0, std::memory_order_relaxed);
	}
	for (int r = 0; r < rn.R; r++) {
		for (int k = 0; k < ROUTER_STATS_SIZE; k++)
			rn.receivers[r]->ctr[k] = 0;
		rn.inbox_pushed[r] = 0;
	}
	for (int p = 0; p < rn.n_nodes; p++) {
		rn.seq_sent[p] = 0;
		rn.n_data_recv[p] = 0;
		rn.end_seq[p] = -1;
		rn.end_sent[p] = 0;
	}
	rn.input_closed.store(0, std::memory_order_relaxed);
	rn.quiescent = false;
	rn.phase = ROUTER_OPEN;
	rn.flush_d = 0;
	rn.flush_built = false;
	rn.flush_blk = ROUTER_NONE;
	rn.flush_k = 0;
	rn.flush_install = false;
	rn.flush_len = 0;
	rn.flush_off = 0;
	rn.round += 1;
	MPI_Barrier(rn.comm);
}

}
#endif
