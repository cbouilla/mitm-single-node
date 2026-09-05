#ifndef MITM_DIRECT_COMMON
#define MITM_DIRECT_COMMON

#include <mpi.h>
#include <vector>
#include <memory>
#include <cmath>
#include <cassert>

#include "tools.hpp"
#include "parameters.hpp"
#include "comm.hpp"

/*
 * The direct scheme's data (PROTOCOL.md §8): its parameters and round header, its counters, and the
 * Scheme that hands them to the engine core.  The two threads are direct_producer.hpp and
 * direct_dict.hpp; the Scheme's functions, the problem wrappers and the entry points are direct.hpp.
 */

namespace mitm::direct {

/* the two phases of a direct round (PROTOCOL.md §8): f fills the dictionary, then g probes it */
enum phase {FILL, PROBE};


/********************************* parameters ********************************/

/*
 * What the direct scheme needs beyond the core's Parameters: the domain, how much of it a round
 * inserts, how many rounds that makes, and the slot layout.  Data only.
 */
struct Params : Parameters {
	int n;                                 /* domain bits: a preimage is n-bit */
	int m;                                 /* range bits: an image, the routing key, is m-bit */
	u64 domain;                            /* 2^n */
	u64 per_round;                         /* entries a round inserts: fill * w, at most the domain */
	u64 n_rounds;                          /* ceil(domain / per_round), capped by --nrounds */
	bool capped;                           /* --nrounds cut the round count: the search is not exhaustive */
	int check_bits;                        /* bits of the key a slot keeps beside the preimage: 63 - n */

	/* n, m are the wrapped problem's.  No producer-per-dict rule: a match is resolved where it is found */
	Params(const Options &o, u64 nbytes_memory, int n, int m) : Parameters(o, nbytes_memory, false), n(n), m(m)
	{
		if (n > 63)
			errx(1, "direct: a %d-bit preimage leaves an 8-byte slot no room for its occupancy bit (n <= 63)", n);
		if (m > 64)
			errx(1, "direct: %d-bit images do not fit a 64-bit key", m);
		if (not (fill > 0 && fill <= 0.9))
			errx(1, "--fill %.2f: the dictionary fill ratio must be in (0, 0.9]", fill);
		domain = 1ull << n;
		per_round = fill * (double) w;
		if (per_round == 0)
			errx(1, "RAM budget too small: %.2f * %" PRIu64 " slots holds no entry", fill, w);
		if (per_round > domain)
			per_round = domain;
		n_rounds = (domain + per_round - 1) / per_round;
		capped = (max_versions < n_rounds);
		if (capped)
			n_rounds = max_versions;
		check_bits = 63 - n;

		/* the FILL phase is the short one: pace the reports on it (PROTOCOL.md §2.2) */
		report_points = per_round / ((u64) reports_per_round * n_nodes);
		if (report_points < 1)
			report_points = 1;
	}
};


/*********************************** header **********************************/

/* the round header (PROTOCOL.md §8): sequenced by rank 0, (round, FILL), (round, PROBE), (round + 1, FILL), ... */
struct Header {
	u64 round = 0xffffffffffffffffull;     /* the direct round; "before the first" is what wraps to 0 */
	u64 phase = PROBE;                     /* FILL or PROBE; PROBE here so that the first next is (0, FILL) */
};


/********************************** counters **********************************/

/*
 * Every tally of the scheme.  One u64[N_COUNTERS] per thread (ThreadContext::ctr) is also the wire
 * layout of a progress report (deltas) and of the end-of-round MPI_Reduce (totals): PROTOCOL.md §2.2.
 */
enum counter {
	/* producers */
	N_EVAL = 0,             /* evaluations of f (FILL) or g (PROBE) */
	N_POINTS,               /* points shipped: one per evaluation */
	/* dict threads */
	N_INSERT,               /* entries inserted (FILL) */
	N_PROBE,                /* probes retired (PROBE) */
	N_STEPS,                /* slots visited by inserts and probes: the cost of linear probing */
	N_MATCH,                /* slots whose check bits matched a probe */
	BAD_MATCH,              /* ... but whose preimage does not map to the key: a check-bit false positive */
	N_COLLISIONS,           /* ... and whose preimage does: f(x) == g(y), golden or not */
	/* comm thread */
	STALL_OUT,              /* turns that ended holding points it could not place (PROTOCOL.md §5) */
	STALL_IN,               /* turns that ended with a received buffer or a bucket waiting on a full dict ring */
	N_COUNTERS
};


/******************************** statistics *********************************/

/* nothing beyond the counters: the three types the core expects, empty */
struct ThreadStats {
	ThreadStats(int) {}
};

struct RoundStats {
	void fold(const RoundStats &) {}

	template <class Ctx>
	void collect(std::vector<std::unique_ptr<Ctx>> &) {}

	void reduce(MPI_Comm, int) {}
};

struct Shared {
	Shared(const Params &) {}
};


class DirectDict;   /* direct_dict.hpp; complete in run(), where the one SharedContext is built and destroyed */

/*
 * What the direct scheme hands to the engine core: its types, its constants and the functions the core
 * calls (PROTOCOL.md §7, §8).  Declared here so that the two thread headers can name the contexts built
 * on it; the functions are defined in direct.hpp, after them.
 */
struct Scheme {
	using Params = direct::Params;
	using Header = direct::Header;
	using Dict = DirectDict;
	using ThreadStats = direct::ThreadStats;
	using RoundStats = direct::RoundStats;
	using Shared = direct::Shared;
	static constexpr int N_COUNTERS = direct::N_COUNTERS;
	static constexpr int PACING = N_POINTS;            /* the counter whose volume paces the reports (§2.2) */
	static constexpr bool LOSSLESS = true;             /* a lost point is a missing entry or probe: drop nothing (§5) */
	static constexpr int STALL_OUT = direct::STALL_OUT;   /* the comm thread's tallies (§5) */
	static constexpr int STALL_IN = direct::STALL_IN;

	/* rank 0 sequences the phases; `stop` comes in as the controller's and goes out raised after the last one */
	template <class Wrapper>
	static void next_header(const Params &params, const Wrapper &wrapper, PRNG &prng, Header &h, u64 &stop);

	/* dict thread d builds its shard, once pinned: first touch (§4.1) */
	static void build_dict(SharedContext<Scheme> &shared, const Params &params, int d);

	/* the two worker threads' rounds: direct_producer.hpp and direct_dict.hpp */
	template <class Wrapper>
	static void producer_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
	                            SharedContext<Scheme> &shared, int index);
	template <class Wrapper>
	static void dict_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
	                        SharedContext<Scheme> &shared, int index);

	/* dict thread d, after the epilogue barrier, given the round's header: the shard is emptied after a PROBE */
	static void after_round(SharedContext<Scheme> &shared, const Params &params, int d, const Header &h);

	/* the controller never closes a round on volume: a phase ends when every producer is exhausted (§4.5) */
	static bool round_complete(const Params &params, const u64 reported[]);
	static void banner(const Params &params, u64 seed);
	static void display(const Params &params, const u64 reported[], double delta, u64 nround);
	static void round_report(const Params &params, const u64 r[], const u64 total[], const RoundStats &round,
	                         const RoundStats &all, double delta, u64 nround);
	static void done(const Params &params, u64 nround, bool found, double seconds);
};

}
#endif
