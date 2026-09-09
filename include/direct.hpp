#ifndef MITM_DIRECT
#define MITM_DIRECT

#include <mpi.h>
#include <omp.h>
#include <err.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <vector>

#include "tools.hpp"
#include "problem.hpp"
#include "parameters.hpp"
#include "router/router.hpp"

/*
 * The direct engine: the exhaustive meet-in-the-middle over a distributed dictionary, on the Router.
 * ceil(2^n / w') rounds of two phases each, w' = fill * w entries: FILL inserts f on the round's chunk of the
 * domain, PROBE probes g on the whole domain.  One Router round per phase, so the Router's own end of round
 * is the barrier the algorithm needs -- every entry of the round is in its shard before the first probe of
 * the round is delivered -- and nothing is dropped, so "no solution" is a proof.
 */

namespace mitm::direct {

using fmt::print;       /* the reports format with {fmt}, unqualified */

/* the two phases of a direct round: f fills the dictionary, then g probes it */
enum phase {FILL, PROBE};

/* the tag of every Router message on the communicator; nothing else on it carries this tag */
static constexpr int ROUTER_TAG = 1;


/********************************* parameters ********************************/

/*
 * Everything a run needs beyond the Options, derived once before the team and const from then on: the MPI
 * topology, the team's shape, the dictionary's slots, the domain and how many rounds cover it.  Data only.
 */
struct Params : Options {
	int rank;                              /* this node: one rank per node */
	int n_nodes;                           /* MPI ranks */
	int S;                                 /* producers per node: the Router's senders */
	int R;                                 /* dict threads per node: the Router's receivers */
	int n_threads;                         /* 1 + R + S: the OpenMP team */
	int n_producers;                       /* producers over all nodes */
	int n_dicts;                           /* dict threads over all nodes */
	u64 w;                                 /* slots in the whole, distributed dictionary */
	u64 w_shard;                           /* slots per shard */
	int n;                                 /* domain bits: a preimage is n-bit */
	u64 domain;                            /* 2^n */
	u64 per_round;                         /* entries a round inserts: fill * w, at most the domain */
	u64 n_rounds_full;                     /* rounds an exhaustive search needs: ceil(domain / per_round) */
	u64 n_rounds;                          /* rounds this run will do: n_rounds_full, capped by --nrounds */
	bool capped;                           /* --nrounds cut the round count: the search is not exhaustive */

	Params(const Options &o, u64 nbytes_memory, int n, int m) : Options(o), n(n)
	{
		MPI_Comm_rank(mpi_comm, &rank);
		MPI_Comm_size(mpi_comm, &n_nodes);
		router.verbose = router.verbose && verbose;   /* the Router prints on rank 0 by itself */
		verbose = verbose && (rank == 0);
		if (n > 63)
			errx(1, "direct: a %d-bit preimage leaves an 8-byte slot no room for its occupancy bit (n <= 63)", n);
		if (m > 64)
			errx(1, "direct: %d-bit images do not fit a 64-bit key", m);
		if (dicts_per_node < 1)
			errx(1, "--dicts-per-node %d: at least one dictionary shard per node", dicts_per_node);
		R = dicts_per_node;
		S = producers_per_node;
		int n_cpu = omp_get_num_procs();   /* the CPUs the rank may use: what the team fills when W == 0 */
		if (S == 0)
			S = n_cpu - 1 - R;
		if (S < 1)
			errx(1, "direct: no CPU left for a producer (%d available, 1 service thread, %d dict threads)",
			     n_cpu, R);
		producers_per_node = S;
		n_threads = 1 + R + S;
		n_producers = n_nodes * S;
		n_dicts = n_nodes * R;

		if (nbytes_memory == 0)
			errx(1, "the RAM budget per node must be given (and nonzero)");
		w_shard = (nbytes_memory * n_nodes / sizeof(u64)) / n_dicts;
		w = w_shard * n_dicts;
		if (w_shard == 0)
			errx(1, "RAM budget too small: %" PRIu64 " bytes/node cannot hold one slot per shard", nbytes_memory);
		if (not (fill > 0 && fill <= 0.9))
			errx(1, "--fill %.2f: the dictionary fill ratio must be in (0, 0.9]", fill);
		domain = 1ull << n;
		per_round = fill * (double) w;
		if (per_round == 0)
			errx(1, "RAM budget too small: %.2f * %" PRIu64 " slots holds no entry", fill, w);
		if (per_round > domain)
			per_round = domain;
		n_rounds_full = (domain + per_round - 1) / per_round;
		n_rounds = n_rounds_full;
		capped = (max_versions < n_rounds);
		if (capped)
			n_rounds = max_versions;
	}
};


/********************************** counters **********************************/

/* every tally of the scheme; one u64[N_COUNTERS] per thread, written by its owner, summed by thread 0 */
enum counter {
	N_EVAL = 0,             /* producers: evaluations of f (FILL) or g (PROBE), one point pushed each */
	N_INSERT,               /* dict threads: entries inserted (FILL) */
	N_PROBE,                /* probes retired (PROBE) */
	N_STEPS,                /* slots visited by inserts and probes: the cost of linear probing */
	N_MATCH,                /* slots whose check bits matched a probe */
	BAD_MATCH,              /* ... but whose preimage does not map to the key: a check-bit false positive */
	N_COLLISIONS,           /* ... and whose preimage does: f(x) == g(y), golden or not */
	N_COUNTERS
};

/* the round record a node contributes to the epilogue's MPI_Allgather: counters, Router stats, its golden pair */
enum record {
	REC_ROUTER = N_COUNTERS,                        /* ROUTER_STATS_SIZE words: the node's Router_Stats */
	REC_FOUND = N_COUNTERS + ROUTER_STATS_SIZE,     /* 1 if the node holds a golden pair */
	REC_X,                                          /* its x */
	REC_Y,                                          /* its y */
	REC_WORDS
};


/********************************** shared state **********************************/

/* one thread's tallies, on a cache line of its own: written by the owner, read by thread 0 */
struct alignas(64) Tally {
	u64 ctr[N_COUNTERS];                   /* indexed by enum counter */
};

/* what a rank's threads share: the tallies, the golden pair a dict thread found, the epilogue's verdict */
struct Shared {
	std::vector<Tally> tally;              /* per OpenMP thread id */
	std::mutex golden_mtx;                 /* serialises set_golden */
	std::atomic<u32> found{0};             /* 1 once golden[] holds this node's pair */
	u64 golden[2] = {};                    /* x and y of the first golden pair found on this node */
	u64 stop = 0;                          /* no next phase; thread 0 writes it, everyone reads it after the barrier */
	u64 phases = 0;                        /* phases run to their end, thread 0's count */
	bool solved = false;                   /* the search found a pair, on this node or another */
	u64 solution[2] = {};                  /* x and y, the same on every node */

	Shared(int n_threads) : tally(n_threads) {}

	/* a dict thread's golden pair; the first one wins */
	void set_golden(u64 x, u64 y)
	{
		std::lock_guard<std::mutex> lock(golden_mtx);
		if (found.load(std::memory_order_relaxed))
			return;
		golden[0] = x;
		golden[1] = y;
		found.store(1, std::memory_order_release);
	}
};


/******************************** the problem wrappers ********************************/

/*
 * A claw problem, f, g : {0,1}^n -> {0,1}^m.  FILL evaluates f, PROBE evaluates g; a match (x, y) is
 * verified with f(x) and tested with is_good_pair(x, y).
 */
template <class Problem>
class ClawWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	bool choice[2][Problem::vlen];  /* what vfg selects, one vector per phase: all f in FILL, all g in PROBE */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "claw";    /* what the banner calls this search */
	static constexpr const char *funcs = "f, g";   /* the functions it evaluates */

	ClawWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n))
	{
		static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
		              "problem not derived from mitm::AbstractClawProblem");
		assert(n <= 63 && m <= 64);
		for (int k = 0; k < vlen; k++) {
			choice[FILL][k] = true;
			choice[PROBE][k] = false;
		}
	}

	/* vlen evaluations of the phase's function: f while filling, g while probing */
	void veval(int phase, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1)
			y[0] = (phase == FILL) ? pb.f(x[0]) : pb.g(x[0]);
		else
			pb.vfg(x, choice[phase], y);
	}

	/* is the collision f(x) == g(y) the one we want? */
	bool good(u64 x, u64 y) const
	{
		return pb.is_good_pair(x, y);
	}

	/* one value every rank must agree on: the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.g(b & in_mask);
	}
};

/*
 * A collision problem, f : {0,1}^n -> {0,1}^m: both phases evaluate f.  A match is a pair of distinct
 * preimages, which is_good_pair -- symmetric, by contract -- accepts or not.
 */
template <class Problem>
class CollisionWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "collision";   /* what the banner calls this search */
	static constexpr const char *funcs = "f";          /* the function it evaluates */

	CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n))
	{
		static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
		              "problem not derived from mitm::AbstractCollisionProblem");
		assert(n <= 63 && m <= 64);
	}

	void veval(int, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1)
			y[0] = pb.f(x[0]);
		else
			pb.vf(x, y);
	}

	/* is the collision f(x) == f(y) the one we want?  is_good_pair is symmetric, so either order will do */
	bool good(u64 x, u64 y) const
	{
		return x != y && pb.is_good_pair(x, y);
	}

	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.f(b & in_mask);
	}
};


/******************************* dictionary shard *****************************/

/*
 * One shard of the distributed dictionary: linear probing over 8-byte slots, insert-only while the round fills,
 * probe-only while it probes, emptied after that.  A slot is the preimage in its low n bits, the key's check
 * bits above, and bit 63 set; zero is empty.  The key is the hashed image the Router delivered, mixed once more
 * here: the routing consumed the top bits of its low word, and a shard above 2^32 slots cut from the key itself
 * would leave part of its slots unreachable.  Built by its dict thread once it is pinned.
 */
class DirectDict {
public:
	static constexpr u64 OCCUPIED = 1ull << 63;   /* the occupancy bit, above the check bits */

	const u64 n_slots;     /* size of A */
	const int n;           /* preimage bits */
	const u64 xmask;       /* the low n bits: the preimage */
	const u64 cmask;       /* the check bits, 63 - n of them, before their shift */
	vector<u64> A;         /* A[i] == OCCUPIED | check << n | preimage, or 0 */

	DirectDict(u64 n_slots, int n) : n_slots(n_slots), n(n), xmask(make_mask(n)), cmask(make_mask(63 - n))
	{
		assert(n <= 63);
		A.resize(n_slots);     /* zero-filled by the calling thread: the first touch of every page */
	}

	/*
	 * Store (key, x) in the first empty slot of its run, which starts at the top bits of the mixed key, scaled
	 * to the table by one multiplication.  Tallies the insert and the slots it visited.
	 */
	void insert(u64 key, u64 x, u64 *ctr)
	{
		u64 h = murmur64(key);
		u64 slot = OCCUPIED | ((h & cmask) << n) | x;
		u64 s = (u64) (((unsigned __int128) h * n_slots) >> 64);
		for (u64 k = 0; k < n_slots; k++, s = (s + 1 == n_slots) ? 0 : s + 1) {
			if (A[s] != 0)
				continue;
			A[s] = slot;
			ctr[N_INSERT] += 1;
			ctr[N_STEPS] += k;
			return;
		}
		errx(1, "direct: a dictionary shard is full (%" PRIu64 " slots): lower --fill", n_slots);
	}

	/*
	 * One probe (key, y) against the shard: every slot of the key's run carrying its tag -- bit 63 and the
	 * mixed key's low check bits -- is a match, kept only if the preimage really maps to the key, and a true
	 * collision is tested for the golden pair, which goes to `shared`.
	 */
	template <class Wrapper>
	void probe(const Wrapper &wrapper, Shared &shared, u64 *ctr, u64 round, u64 key, u64 y) const
	{
		u64 h = murmur64(key);
		u64 tag = OCCUPIED | ((h & cmask) << n);
		u64 s = (u64) (((unsigned __int128) h * n_slots) >> 64);
		u64 k = 0;
		for (; k < n_slots && A[s] != 0; k++, s = (s + 1 == n_slots) ? 0 : s + 1) {
			if ((A[s] & ~xmask) != tag)
				continue;
			ctr[N_MATCH] += 1;
			u64 x = A[s] & xmask;
			if (murmur64(wrapper.pb.f(x)) != key) {
				ctr[BAD_MATCH] += 1;
				continue;
			}
			ctr[N_COLLISIONS] += 1;
			if (not wrapper.good(x, y))
				continue;
			printf("\nFound golden pair! round=%" PRIu64 " x=%" PRIx64 " y=%" PRIx64 "\n", round, x, y);
			shared.set_golden(x, y);
		}
		ctr[N_PROBE] += 1;
		ctr[N_STEPS] += k;
	}

	/* empty the shard after a PROBE phase; the owner alone today */
	void flush()
	{
		memset(A.data(), 0, n_slots * sizeof(u64));
	}
};


/******************************* the dict thread ******************************/

/*
 * A dict thread's phase: take the points the Router delivers, one at a time, until nothing more can arrive.
 * While the round fills, a point is an entry, (hashed image, preimage), inserted; while it probes, a point is a
 * probe, resolved on the spot: a match costs about what a probe does, so nothing is handed to anybody.  After a
 * PROBE phase the shard is emptied.
 */
template <class Wrapper>
void dict_round(Router_thread &rt, const Wrapper &wrapper, Shared &shared, DirectDict &dict, u64 *ctr,
                u64 round, int phase)
{
	for (;;) {
		u64 key;                       /* the hashed image the producer routed on */
		u64 val;                       /* its preimage */
		if (not Router_Pop(&key, &val, rt)) {
			if (Router_Test_drained(rt))
				break;
			cpu_relax();               /* no CAS in Router_Pop --> no need for backoff */
			continue;
		}
		if (phase == FILL)
			dict.insert(key, val, ctr);
		else
			dict.probe(wrapper, shared, ctr, round, key, val);
	}
	if (phase == PROBE)
		dict.flush();
}


/******************************** the producer thread ********************************/

/*
 * A producer's phase: evaluate the phase's function on its own contiguous piece of the phase's span -- the
 * round's chunk while filling, the whole domain while probing -- vlen inputs at a time, push every image as
 * (murmur64(image), preimage) to the receiver the hash's low word names, then close.  The pieces of the
 * n_producers producers cut the span into equal parts, the same on every node, so the partition needs no
 * message; the hash spares the producer a 64-bit division per point.
 */
template <class Wrapper>
void producer_round(Router_thread &rt, const Wrapper &wrapper, const Params &params, u64 *ctr, u64 round, int phase)
{
	constexpr int vlen = Wrapper::vlen;
	const u64 n_recv = Router_num_recv(rt);
	u64 lo = 0;
	u64 hi = params.domain;
	if (phase == FILL) {
		lo = round * params.per_round;
		hi = std::min(lo + params.per_round, params.domain);
	}
	u64 p = Router_rank(rt);
	u64 n_pieces = Router_num_send(rt);
	u64 span = hi - lo;
	u64 piece = span / n_pieces;
	u64 extra = span % n_pieces;
	u64 my_lo = lo + piece * p + std::min(p, extra);
	u64 my_hi = my_lo + piece + ((p < extra) ? 1 : 0);

	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
	for (u64 base = my_lo; base < my_hi; base += vlen) {
		int valid = (my_hi - base < (u64) vlen) ? (int) (my_hi - base) : vlen;
		for (int k = 0; k < vlen; k++)
			x[k] = (k < valid) ? base + k : base;   /* a lane past the end still gets an input in the domain */
		wrapper.veval(phase, x, y);
		for (int k = 0; k < valid; k++) {
			u64 h = murmur64(y[k]);
			int dest = (int) (((h & 0xffffffffull) * n_recv) >> 32);
			Router_Push(h, x[k], dest, rt);
		}
		ctr[N_EVAL] += valid;
	}
	Router_Close(rt);
}


/******************************** the printing ********************************/

/* the startup report, on rank 0 before the team exists; the Router prints its own layout after it */
template <class Wrapper>
static void banner(const Wrapper &wrapper, const Params &params, u64 seed)
{
	print("Starting direct {} search with {} : {{0,1}}^{} --> {{0, 1}}^{} (vlen={})\n", Wrapper::kind,
	      Wrapper::funcs, wrapper.n, wrapper.m, Wrapper::vlen);
	print("Starting the direct (exhaustive) meet-in-the-middle on the Router with seed={:016x}\n", seed);
	print("MPI: {} node(s) x (1 service + {} dict + {} producer) = {} threads/node\n", params.n_nodes,
	      params.R, params.S, params.n_threads);
	print("MPI: {} dictionary shards, {} producer threads in total\n", params.n_dicts, params.n_producers);
	double work = (double) params.n_rounds * (params.per_round + params.domain);
	std::string hper = human_format(params.per_round);
	print("Dictionary: linear probing, 8-byte slots = {}-bit preimage | {} check bits | occupancy.  "
	      "Fill {:.2f}: {} entries per round\n", params.n, 63 - params.n, params.fill, hper);
	print("RAM per node == {:.1f} MB of dictionary; {} slots in all (2^{:.2f}), {} per shard\n",
	      (double) params.w_shard * params.R * sizeof(u64) / 1e6, human_format(params.w),
	      std::log2((double) params.w), hper);
	print("{} round(s) of two phases: FILL (f on 2^{:.2f} preimages) then PROBE (g on all 2^{}).  "
	      "Total work: {} evaluations = 2^{:.2f}\n", params.n_rounds, std::log2((double) params.per_round),
	      params.n, human_format(work), std::log2(work));
	if (params.capped)
		print("NOTICE: --nrounds caps the search at {} round(s): it is NOT exhaustive\n",
		      params.n_rounds);
	fflush(stdout);
}

/* the live one-line refresh: rank 0's own node, read from its tallies without synchronisation */
static void display(const Params &params, const Shared &shared, const u64 *stats, double delta, u64 round,
                    int phase)
{
	u64 eval = 0;
	u64 retired = 0;
	for (int t = 0; t < params.n_threads; t++) {
		eval += shared.tally[t].ctr[N_EVAL];
		retired += shared.tally[t].ctr[N_INSERT] + shared.tally[t].ctr[N_PROBE];
	}
	u64 span = params.domain;
	if (phase == FILL)
		span = std::min(params.per_round, params.domain - round * params.per_round);
	double completion = (double) eval * params.n_nodes / (double) span;
	print("\rRound {}/{} {}:  {:.1f}s ({:.1f}%, ETA {:.1f}s).  {} #f/s per producer.  "
	      "{} points/s per dict thread.  node-->{}B/s   ", round, params.n_rounds,
	      (phase == FILL) ? "FILL" : "PROBE", delta, 100. * completion,
	      (completion > 0) ? delta * (1 - completion) / completion : 0.,
	      human_format((double) eval / params.S / delta),
	      human_format((double) retired / params.R / delta),
	      human_format((double) stats[ROUTER_BYTES_SENT] / delta));
	fflush(stdout);
}

/* the end-of-phase report: `r` is the exact sum over every thread of every node, `total` the all-time sums */
static void round_report(const Params &params, const u64 r[], const u64 total[], double delta, u64 round,
                         int phase)
{
	print("\n");
	print("Round {} {}.  {:.1f}s.  2^{:.2f} evaluations (total 2^{:.2f}).  {} #f/s per producer.  "
	      "node-->{}B/s.  {} points/s per dict thread\n", round, (phase == FILL) ? "FILL" : "PROBE", delta,
	      std::log2((double) (r[N_EVAL] ? r[N_EVAL] : 1)),
	      std::log2((double) (total[N_EVAL] ? total[N_EVAL] : 1)),
	      human_format((double) r[N_EVAL] / params.n_producers / delta),
	      human_format((double) r[REC_ROUTER + ROUTER_BYTES_SENT] / params.n_nodes / delta),
	      human_format((double) (r[N_INSERT] + r[N_PROBE]) / params.n_dicts / delta));
	if (r[N_INSERT] > 0)
		print("            {} inserted, load {:.2f}/slot, {:.2f} slots visited per insert\n",
		      r[N_INSERT], (double) r[N_INSERT] / params.w, (double) r[N_STEPS] / r[N_INSERT]);
	if (r[N_PROBE] > 0)
		print("            {} probed, {:.2f} slots visited per probe.  {} matches: {} false positives, "
		      "{} collisions (total 2^{:.2f})\n", r[N_PROBE], (double) r[N_STEPS] / r[N_PROBE],
		      r[N_MATCH], r[BAD_MATCH], r[N_COLLISIONS],
		      std::log2((double) (total[N_COLLISIONS] ? total[N_COLLISIONS] : 1)));
	if (r[REC_ROUTER + ROUTER_STALL_OUT] | r[REC_ROUTER + ROUTER_STALL_IN])
		print("            STALLED  {} blocks held back for a full destination / {} receives left "
		      "unposted for want of a free block\n", r[REC_ROUTER + ROUTER_STALL_OUT],
		      r[REC_ROUTER + ROUTER_STALL_IN]);
	if (r[REC_ROUTER + ROUTER_PUSHED] != r[REC_ROUTER + ROUTER_POPPED])
		print("            ACCOUNTING BROKEN: {} points pushed, {} delivered\n",
		      r[REC_ROUTER + ROUTER_PUSHED], r[REC_ROUTER + ROUTER_POPPED]);
	print("\n");
	fflush(stdout);
}

/* the last line: found, or the domain exhausted -- a proof of absence, unless --nrounds cut it short */
static void done(const Params &params, u64 phases, bool found, double seconds)
{
	u64 rounds = (phases + 1) / 2;
	if (found)
		print("Completed in {:.2f}s\n", seconds);
	else if (params.capped)
		print("Gave up after {} of the {} round(s) an exhaustive search needs ({:.2f}s)\n",
		      rounds, params.n_rounds_full, seconds);
	else
		print("No solution: the domain is exhausted after {} round(s) ({:.2f}s)\n", rounds, seconds);
	fflush(stdout);
}


/******************************** the service thread and the epilogue ********************************/

/*
 * Thread 0's phase: turn the Router until the node has sent and received everything, refreshing the live
 * line, then read the node's Router tallies, which Router_Reset clears.
 */
static void service_round(Router_thread &rt, const Params &params, const Shared &shared, u64 *stats, u64 round,
                          int phase, double t0)
{
	double last = t0;
	u64 turns = 0;
	while (not Router_Test_quiescent(rt)) {
		Router_Progress(rt);
		turns += 1;
		if (not params.verbose || (turns & 0xff) != 0)
			continue;
		double now = wtime();
		if (now - last < params.ping_delay)
			continue;
		last = now;
		Router_Stats(stats, rt);
		display(params, shared, stats, now - t0, round, phase);
	}
	Router_Stats(stats, rt);
}

/*
 * Thread 0's end of phase, once the Router is reset: the node's record -- its tallies, its Router stats and
 * the golden pair it may hold -- goes into one Allgather, and every node reads the same verdict out of it.
 * The lowest rank that found a pair provides the answer.
 */
static void epilogue(const Params &params, Shared &shared, const u64 *stats, u64 *records, u64 *total,
                     u64 round, int phase, double delta)
{
	u64 rec[REC_WORDS] = {};
	for (int t = 0; t < params.n_threads; t++)
		for (int c = 0; c < N_COUNTERS; c++)
			rec[c] += shared.tally[t].ctr[c];
	for (int k = 0; k < ROUTER_STATS_SIZE; k++)
		rec[REC_ROUTER + k] = stats[k];
	if (shared.found.load(std::memory_order_acquire)) {
		rec[REC_FOUND] = 1;
		rec[REC_X] = shared.golden[0];
		rec[REC_Y] = shared.golden[1];
	}
	MPI_Allgather(rec, REC_WORDS, MPI_UINT64_T, records, REC_WORDS, MPI_UINT64_T, params.mpi_comm);

	u64 sum[REC_WORDS] = {};
	for (int r = 0; r < params.n_nodes; r++) {
		const u64 *q = records + (size_t) r * REC_WORDS;
		for (int k = 0; k < REC_FOUND; k++)
			sum[k] += q[k];
		if (q[REC_FOUND] && not shared.solved) {
			shared.solved = true;
			shared.solution[0] = q[REC_X];
			shared.solution[1] = q[REC_Y];
		}
	}
	for (int k = 0; k < REC_FOUND; k++)
		total[k] += sum[k];
	shared.phases += 1;
	shared.stop = shared.solved || (phase == PROBE && round + 1 >= params.n_rounds);
	if (params.verbose)
		round_report(params, sum, total, delta, round, phase);
}

/******************************** the engine ********************************/

/*
 * The search itself: one OpenMP team per node -- thread 0 the Router's service thread, the next R its dict
 * threads, the rest its producers -- running the same deterministic sequence of phases, (0, FILL), (0, PROBE),
 * (1, FILL), ..., each one Router round.  Returns x and y of the golden pair, the same on every node, or
 * nothing once the domain is exhausted.
 */
template <class Wrapper>
optional<pair<u64, u64>> run(const Wrapper &wrapper, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	const Params params(opts, nbytes_memory, wrapper.n, wrapper.m);

	/* self-test */
	u64 mine[3];
	mine[0] = prng.rand();
	mine[1] = prng.rand();
	mine[2] = wrapper.self_test(mine[0], mine[1]);
	u64 theirs[3] = {mine[0], mine[1], mine[2]};
	MPI_Bcast(theirs, 3, MPI_UINT64_T, 0, params.mpi_comm);
	if (theirs[0] != mine[0] || theirs[1] != mine[1] || theirs[2] != mine[2])
		errx(1, "direct: the ranks do not hold the same problem instance (self-test mismatch)");

	if (params.verbose)
		banner(wrapper, params, prng.seed);

	Shared shared(params.n_threads);
	double t_start = wtime();

	#pragma omp parallel num_threads(params.n_threads)
	{
		int tid = omp_get_thread_num();
		int role = ROUTER_SENDER;
		if (tid == 0)
			role = ROUTER_SERVICE;
		else if (tid <= params.R)
			role = ROUTER_RECEIVER;
		Router_thread rt = Router_Init(role, ROUTER_GROUP_AUTO, params.mpi_comm, ROUTER_TAG, &params.router);
		u64 *ctr = shared.tally[tid].ctr;

		/* a dict thread's shard, zero-filled here: its own first touch.  Empty on the other roles */
		DirectDict dict(role == ROUTER_RECEIVER ? params.w_shard : 0, params.n);

		std::vector<u64> records;              /* thread 0: the nodes' records, from the Allgather */
		std::vector<u64> total;                /* thread 0: the all-time sums */
		std::vector<u64> stats;                /* thread 0: the node's Router tallies for the phase */
		if (role == ROUTER_SERVICE) {
			if (Router_num_recv(rt) != params.n_dicts)
				errx(1, "direct: the Router has %d receivers, the dictionary %d shards",
				     Router_num_recv(rt), params.n_dicts);
			records.assign((size_t) params.n_nodes * REC_WORDS, 0);
			total.assign(REC_WORDS, 0);
			stats.assign(ROUTER_STATS_SIZE, 0);
		}

		for (u64 step = 0;; step++) {
			u64 round = step / 2;
			int phase = (step % 2 == 0) ? FILL : PROBE;
			for (int c = 0; c < N_COUNTERS; c++)
				ctr[c] = 0;
			double t0 = wtime();

			if (role == ROUTER_SERVICE)
				service_round(rt, params, shared, stats.data(), round, phase, t0);
			else if (role == ROUTER_RECEIVER)
				dict_round(rt, wrapper, shared, dict, ctr, round, phase);
			else
				producer_round(rt, wrapper, params, ctr, round, phase);

			Router_Reset(rt);      /* its own team barriers and MPI_Barrier: every tally is written by now */

			if (role == ROUTER_SERVICE)
				epilogue(params, shared, stats.data(), records.data(), total.data(), round, phase, wtime() - t0);

			#pragma omp barrier    /* the verdict, and every thread's next tally, are thread 0's to publish */
			if (shared.stop)
				break;
		}
	}

	if (params.verbose)
		done(params, shared.phases, shared.solved, wtime() - t_start);
	if (not shared.solved)
		return nullopt;
	return optional(pair(shared.solution[0], shared.solution[1]));
}


/******************************** the entry points ********************************/

/* find x0 != x1 with f(x0) == f(x1) and is_good_pair(x0, x1), or prove there is none */
template <class Problem>
optional<pair<u64, u64>> collision_search(const Problem &pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	CollisionWrapper<Problem> wrapper(pb);
	auto collision = run(wrapper, nbytes_memory, opts, prng);
	if (not collision)
		return nullopt;

	auto [x0, x1] = *collision;
	assert(x0 != x1);
	assert(pb.f(x0) == pb.f(x1));
	assert(pb.is_good_pair(x0, x1));
	return optional(pair(x0, x1));
}

/* find x0, x1 with f(x0) == g(x1) and is_good_pair(x0, x1), or prove there is none */
template <class Problem>
optional<pair<u64, u64>> claw_search(const Problem &pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	ClawWrapper<Problem> wrapper(pb);
	auto claw = run(wrapper, nbytes_memory, opts, prng);
	if (not claw)
		return nullopt;

	auto [x0, x1] = *claw;
	assert(pb.f(x0) == pb.g(x1));
	assert(pb.is_good_pair(x0, x1));
	return optional(pair(x0, x1));
}

}
#endif
