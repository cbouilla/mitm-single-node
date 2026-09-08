#ifndef MITM_ROUTER_COMMON
#define MITM_ROUTER_COMMON

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
#include <unistd.h>
#include <sys/syscall.h>

#include "../tools.hpp"
#include "../parameters.hpp"
#include "topology.hpp"
#include "placement.hpp"
#include "ring.hpp"

namespace mitm {

/****************************** the interface's constants ******************************/

/* router_role and ROUTER_GROUP_AUTO are defined in placement.hpp, which both this and RouterPlacement need */

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
	bool pin = true;                    /* pin the threads and form the groups, or leave the CPUs to the caller */
	int cache_level = 0;                /* the cache level a group sits in; 0 == the lowest shared by several cores */
	int group_size = 16;                /* cores per group at most: the AUTO knob that splits a cache domain */
	bool verbose = true;                /* rank 0 prints the sizes once connected */
};

/* Router_Init's `opts` for the defaults */
static constexpr const Router_Opts *ROUTER_DEFAULT_OPTS = NULL;


/****************************** internal constants and wire format ******************************/

static constexpr u32 ROUTER_NONE = 0xffffffffu;   /* "no block": a bare destination, an empty list, an idle slot */
static constexpr size_t ROUTER_HDR_BYTES = 64;    /* a block's header slot: the message header, one cache line */
static constexpr size_t ROUTER_DEST_WORDS = 8;    /* u64 words per destination: its two words fill one cache line */
static constexpr u32 ROUTER_MAX_L = 4096;         /* lines per block at most */
static constexpr u32 ROUTER_BATCH = 32;           /* blocks off the free ring at a time: a sender's cache, the stash */
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

/*
 * Copy into memory the writer never reads back, without fetching the destination for ownership first: worth 2.1x
 * the point rate of memcpy on a 256-core EPYC, where staging a line pulled its block in from DRAM, half of it
 * from the other socket.  The stores are visible to another thread once this returns.
 * `dst` must be 64-byte aligned and `bytes` a multiple of 64, which a line slot of a block always is.
 */
static inline void router_stream_copy(void *dst, const void *src, size_t bytes)
{
#if defined(__AVX512F__)
	__m512i *d = (__m512i *) dst;
	const __m512i *s = (const __m512i *) src;
	for (size_t i = 0; i < bytes / sizeof(*d); i++)
		_mm512_stream_si512(d + i, _mm512_loadu_si512(s + i));
	_mm_sfence();                       /* nothing else orders a non-temporal store before a release store */
#elif defined(__AVX2__)
	__m256i *d = (__m256i *) dst;
	const __m256i *s = (const __m256i *) src;
	for (size_t i = 0; i < bytes / sizeof(*d); i++)
		_mm256_stream_si256(d + i, _mm256_loadu_si256(s + i));
	_mm_sfence();                       /* nothing else orders a non-temporal store before a release store */
#else
	memcpy(dst, src, bytes);            /* no non-temporal store here: the destination is fetched for ownership */
#endif
}

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
 * A sender also holds a cache of free blocks, zeroed in advance: the block it installs when it seals comes out
 * of it, so the install the other senders of that destination wait on is one load and one store, and the free
 * ring sees the sender once in ROUTER_BATCH seals, after the seal, when nobody waits on it.
 * A receiver holds whole blocks: its inbox names them, Router_Grab hands it the next one to read in place, and
 * Router_Release pushes it onto the node's free ring.
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
	const int group;                     /* a sender's or receiver's group; -1 for the service */
	const int domain;                    /* the cache domain it was pinned in; -1 when not pinned */
	const int cpu;                       /* the CPU the kernel reports after pinning; momentary when not pinned */
	const int numa_node;                 /* the NUMA node the kernel reports after pinning */
	Router_node &node;                   /* the node it belongs to */
	/* a sender's cache of free blocks, zeroed already: what it installs at a seal, refilled off the free ring */
	u32 cache[ROUTER_BATCH];             /* the blocks, cache_n of them, the next to install last */
	u32 cache_n = 0;                     /* how many */
	/* a sender's flag, on a line of its own: the service reads it every turn */
	alignas(64) std::atomic<u32> closed; /* release-stored by Router_Close, acquired by the service */
	/* a receiver's */
	RouterRing inbox;                    /* from the service: (block, count), blocks it now holds */
	u32 cur_blk = ROUTER_NONE;           /* the block out: grabbed, not yet released; NONE when none */
	const Point *cur_pts = NULL;         /* its points, in block memory */
	u32 cur_count = 0;                   /* how many */
	u32 cur_off = 0;                     /* Router_Pop's cursor into it */

	Router_thread(int role, int index, int group, int domain, int cpu, int numa_node, Router_node &rn);
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
	std::vector<int> colors;             /* per thread, its group color, or ROUTER_GROUP_AUTO */
	std::vector<cpu_set_t> masks;        /* per thread, its own affinity mask: connect() unions them */
	RouterPlacement plan;                /* the thread layout: the groups, the domains, the pinning */
	int misplaced = 0;                   /* threads the kernel did not pin as planned: Router_Init aborts if any */
	std::vector<Router_thread *> senders;                     /* the S senders' objects (on their stacks), by index */
	std::vector<Router_thread *> receivers;                   /* the R receivers' objects (on their stacks), by index */
	/* a destination's line, alone on a cache line.  Word 0: (lines reserved << 32) | block id, the counter in the
	 * high word so that a runaway counter carries out of the word, not into the block id.  Word 1: the points in
	 * its closing buffer, a closer's fetch_add taking its offset. */
	std::atomic<u64> *dest = NULL;       /* F * ROUTER_DEST_WORDS */
	char *pool = NULL;                   /* n_blocks blocks of block_bytes: staging, send and receive memory alike */
	std::atomic<u8> *n_valid = NULL;     /* n_blocks * nv_stride: byte k is 0 until line k is written, then 1 */
	Point *partial = NULL;               /* F closing buffers of partial_cap points each, filled at Router_Close */
	std::atomic<u64> *blk_link = NULL;   /* per block: its link in the sealed stack, its destination above, or parked */

	/* the free ring: block ids in cells of one word, (sequence << 32) | block.  A pusher's ticket comes off free_in
	 * by fetch_add and never waits for more than the popper of the cell's previous element, the ring being never
	 * full (cells >= n_blocks); a popper claims a ready prefix at free_out by one CAS, or fails at once. */
	std::atomic<u64> *free_cell = NULL;  /* free_mask + 1 cells; cell i starts as (i << 32) | ROUTER_NONE */
	u32 free_mask = 0;                   /* cells - 1, cells the power of two >= n_blocks */
	alignas(64) std::atomic<u64> free_in{0};   /* push tickets issued */
	u64 free_in_pad[7] = {};             /* the rest of that line */
	alignas(64) std::atomic<u64> free_out{0};  /* elements claimed */
	u64 free_out_pad[7] = {};            /* the rest of that line */

	/* the sealed stack, its word alone on a cache line: the sealers push, the service takes it whole */
	alignas(64) std::atomic<u32> sealed_top{ROUTER_NONE};
	u32 sealed_pad[15] = {};             /* the rest of that line */

	/* the service thread's own */
	u64 ctr[ROUTER_STATS_SIZE] = {};     /* its share of the tallies: the network, the drops it makes, the turns */
	std::vector<u32> free_list;          /* its stash: its own frees, for its own needs; a batch to or from the ring */
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

	/* connection: Router_Init; connect and banner on thread 0, touch_pool on every thread */
	void connect();
	void touch_pool(int tid, int n_threads);
	void banner() const;

	/* the free ring: any thread.  Router_Release pushes the block a receiver read through, a sender's cache takes a
	 * batch, the stash trades them */
	void free_push(const u32 *blks, u32 n);
	u32 free_pop_many(u32 k, u32 *out);

	/* the sender path: Router_Push, by the push that fills a private line; refill also from Router_Init */
	void stage_line(Router_thread &s, int d, const Point *line);
	void refill(Router_thread &s);
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

}
#endif
