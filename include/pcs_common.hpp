#ifndef MITM_PCS_COMMON
#define MITM_PCS_COMMON

#include <mpi.h>
#include <mutex>
#include <vector>
#include <memory>
#include <cmath>
#include <cassert>
#include <strings.h>          // ffsll, for the HyperLogLog

#include "tools.hpp"
#include "parameters.hpp"
#include "comm.hpp"

/*
 * The PCS scheme's data (PROTOCOL.md §7): its parameters and round header, its counters, the collision
 * queue a dict thread feeds and the HyperLogLog, and the Scheme that hands all of it to the engine
 * core.  The two threads are walker.hpp and inserter.hpp; the Scheme's functions, the problem wrappers
 * and the entry points are pcs.hpp.
 */

namespace mitm::pcs {

/********************************* parameters ********************************/

/* What PCS needs beyond the core's Parameters: the DP layout and the difficulty.  Data only. */
struct Params : Parameters {
	int jbits;                             /* bits of chain index stored with each DP */
	int lenbits;                           /* width of the length field beside it */
	u64 len_sat;                           /* its largest value == "at least this, exact length unknown" */
	bool theta_auto;                       /* theta was chosen from alpha */
	double auto_theta;                     /* what alpha chooses, whether or not it was used */
	u64 threshold;                         /* x is a DP iff x <= threshold */
	u64 dp_max_it;                         /* how many iterations to find a DP */
	u64 points_per_version;                /* #DP per version of the function */

	/* n, m are the wrapped problem's.  Every dict thread needs a producer: its collision queue (§3.2) */
	Params(const Options &o, u64 nbytes_memory, int n, int m) : Parameters(o, nbytes_memory, true)
	{
		jbits = std::log2(10 * w) + 8;
		if (jbits > 56)
			errx(1, "dictionary too large: a %d-bit chain index leaves a slot no room for a length", jbits);
		lenbits = dp_lenbits ? dp_lenbits : 64 - jbits;
		if (lenbits < 1 || jbits + lenbits > 64)
			errx(1, "--dp-len-bits %d does not fit beside a %d-bit chain index", lenbits, jbits);
		len_sat = make_mask(lenbits);

		auto_theta = alpha * std::sqrt(std::ldexp((double) w, -n));
		theta_auto = (theta < 0);
		if (theta_auto)
			theta = std::min(auto_theta, 1.0);
		threshold = pow(2, m) * theta;
		dp_max_it = 20 / theta;
		points_per_version = beta * w;

		/* report by volume too: on a timer alone a short round overshoots beta*w (PROTOCOL.md §2.2) */
		report_points = points_per_version / ((u64) reports_per_round * n_nodes);
		if (report_points < 1)
			report_points = 1;
	}
};


/*********************************** header **********************************/

/* the round header (PROTOCOL.md §4.2): drawn at random by rank 0 */
struct Header {
	u64 i = 0;                      /* mixing function version */
	u64 root_seed = 0;              /* chain j starts at root_seed + j * multiplier */
};


/********************************** counters **********************************/

/*
 * Every tally of the scheme.  One u64[N_COUNTERS] per thread (ThreadContext::ctr) is also the wire
 * layout of a progress report (deltas) and of the end-of-round MPI_Reduce (totals): PROTOCOL.md §2.2.
 * None of this is part of the attack: it is what tells us whether the parameters are any good.
 */
enum counter {
	/* producers */
	N_EVAL = 0,             /* evaluations of the mixing function, by walking or by resolving */
	N_DP,                   /* distinguished points found */
	N_POINTS_TRAILS,        /* sum of the lengths of the trails that reached a DP */
	N_COLLISIONS,           /* collisions located */
	COLLIDING_LEN_MIN,      /* sum of the shorter length of each colliding pair */
	COLLIDING_LEN_MAX,      /* ... and of the longer one */
	N_MEASURE,              /* trails re-walked because their length had saturated */
	BAD_DP,                 /* trail gave up before a distinguished point, walking or re-walked to be measured */
	BAD_COLLISION,          /* the two trails "collide" on the same value */
	BAD_WALK_ROBINHOOD,     /* one trail is a suffix of the other */
	BAD_WALK_NONCOLLIDING,  /* dictionary false positive: the trails never meet */
	DROP_PRODUCERQ,         /* DP dropped: the queue to the comm thread was full */
	/* dict threads */
	N_PROBE,                /* dictionary probes retired */
	BAD_PROBE,              /* dictionary slot was empty, or held a different key */
	DROP_COLL,              /* candidate dropped: the collision queue was full */
	/* comm thread */
	DROP_OUT,               /* DP dropped: the outgoing MPI buffer was still in flight */
	DROP_DICTQ,             /* DP dropped: a local dict thread's queue was full */
	N_COUNTERS
};


/****************************** collision queue *******************************/

struct CollisionCandidate {
	u64 i;                  /* mixing function version; a producer asserts it is the round's */
	u64 seed0;              /* the incoming point: its chain index, ... */
	u64 len0_maybe;         /* ... its trail length; 0 == it saturated on the wire, the producer re-walks it */
	u64 end;                /* ... and its endpoint, in full: what a re-walked trail must reach */
	u64 seed1;              /* the point that was in the slot: its chain index, ... */
	u64 len1_maybe;         /* ... and its length; 0 == saturated in the dictionary, the producer re-walks it */
};

/*
 * Bounded, mutex-protected queue of dictionary hits: ONE per dict thread, pushed by that dict thread and
 * popped by the producers of its thread group.  Candidates cross it in runs, never one at a time, so
 * that the lock and the `count` line are touched once per run.  Full == the run is truncated and the
 * dict thread tallies the rest (DROP_COLL).  PROTOCOL.md §3.2.
 */
class CollisionQueue {
	std::mutex mtx;
	std::vector<CollisionCandidate> buf;
	size_t head = 0;                /* next to pop */
	size_t tail = 0;                /* next to push */
	std::atomic<size_t> count;      /* probed without the lock on the hot path */

public:
	CollisionQueue(size_t capacity) : buf(capacity < 1 ? 1 : capacity), count(0) {}

	/* Relaxed on purpose: a stale answer costs a producer one wasted pop() at worst. */
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

	/* consumer side.  Several producers of the group may be here at once; the lock is what orders them. */
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


/******************************** HyperLogLog *********************************/

static constexpr int HLL_REGISTERS = 0x10000;   /* one per value of the top 16 bits of a pair's hash */

/*
 * How many DISTINCT collisions a round found.  One HyperLogLog per producer, plain `u8`, written by its
 * owner alone and merged by the comm thread once every producer is quiescent: PROTOCOL.md §3.4.  None of
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

/* a thread's statistics beyond its counters: a producer's HyperLogLog of the round's collisions (§3.4) */
struct ThreadStats {
	std::vector<u8> hll;            /* HLL_REGISTERS plain registers for a producer; empty for the other roles */

	ThreadStats(int role)
	{
		if (role == PRODUCER)
			hll.assign(HLL_REGISTERS, 0);
	}
};

/*
 * A round's statistics beyond its counters: one HyperLogLog -- the producers' merged by the comm thread
 * once every one of them is quiescent, the nodes' MAX-reduced to rank 0, the rounds' folded into the
 * controller's all-time registers (PROTOCOL.md §3.4, §4.6).
 */
struct RoundStats {
	u8 hll[HLL_REGISTERS];

	RoundStats()
	{
		for (int k = 0; k < HLL_REGISTERS; k++)
			hll[k] = 0;
	}

	/* register-wise max, which is why the merge is exact whatever the order */
	void fold(const RoundStats &o)
	{
		hll_merge(hll, o.hll);
	}

	/* merge every producer's registers and zero them for the next round.  Comm thread, workers QUIESCENT */
	template <class Ctx>
	void collect(std::vector<std::unique_ptr<Ctx>> &ctx)
	{
		for (size_t t = 0; t < ctx.size(); t++) {
			if (ctx[t]->role != PRODUCER)
				continue;
			hll_merge(hll, ctx[t]->scheme.hll.data());
			for (int k = 0; k < HLL_REGISTERS; k++)
				ctx[t]->scheme.hll[k] = 0;
		}
	}

	/* MAX to rank 0.  Blocks; every rank reaches it after its drain (§4.6) */
	void reduce(MPI_Comm comm, int rank)
	{
		if (rank == 0)
			MPI_Reduce(MPI_IN_PLACE, hll, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, comm);
		else
			MPI_Reduce(hll, NULL, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, comm);
	}
};


/* the scheme's shared state: the collision queues, one per dict thread, built by it once pinned (§3.2) */
struct Shared {
	std::vector<std::unique_ptr<CollisionQueue>> coll_q;   /* queue r: dict thread r pushes, group r pops */

	Shared(const Params &params) : coll_q(params.dicts_per_node) {}
};


class PcsDict;   /* inserter.hpp; complete in run(), where the one SharedContext is built and destroyed */

/*
 * What PCS hands to the engine core: its types, its constants and the functions the core calls
 * (PROTOCOL.md §7).  Declared here so that walker.hpp and inserter.hpp can name the contexts built on
 * it; the functions are defined in pcs.hpp, after them.
 */
struct Scheme {
	using Params = pcs::Params;
	using Header = pcs::Header;
	using Dict = PcsDict;
	using ThreadStats = pcs::ThreadStats;
	using RoundStats = pcs::RoundStats;
	using Shared = pcs::Shared;
	static constexpr int N_COUNTERS = pcs::N_COUNTERS;
	static constexpr int PACING = N_DP;                /* the counter whose volume paces the reports (§2.2) */
	static constexpr int DROP_OUT = pcs::DROP_OUT;     /* the comm thread's two drop tallies (§5) */
	static constexpr int DROP_DICTQ = pcs::DROP_DICTQ;
	static constexpr bool LOSSLESS = false;            /* a point that finds a full channel is dropped and tallied */

	/* rank 0 draws the next round's header; `stop` comes in as the controller's and goes out final (§4.2) */
	template <class Wrapper>
	static void next_header(const Params &params, const Wrapper &wrapper, PRNG &prng, Header &h, u64 &stop);

	/* dict thread d builds its shard and its collision queue, once pinned: first touch (§4.1) */
	static void build_dict(SharedContext<Scheme> &shared, const Params &params, int d);

	/* the two worker threads' rounds */
	template <class Wrapper>
	static void producer_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
	                            SharedContext<Scheme> &shared, int index);
	template <class Wrapper>
	static void dict_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
	                        SharedContext<Scheme> &shared, int index);

	/* dict thread d, after the epilogue barrier: its shard is emptied for the next round (§4.6) */
	static void after_round(SharedContext<Scheme> &shared, const Params &params, int d);

	/* the controller's one decision, and all of its printing */
	static bool round_complete(const Params &params, const u64 reported[]);
	static void banner(const Params &params, u64 seed);
	static void display(const Params &params, const u64 reported[], double delta, u64 nround);
	static void round_report(const Params &params, const u64 r[], const u64 total[], const RoundStats &round,
	                         const RoundStats &all, double delta, u64 nround);
	static void done(const Params &params, u64 nround, bool found, double seconds);
};

}
#endif
