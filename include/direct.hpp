#ifndef MITM_DIRECT
#define MITM_DIRECT

#include <cassert>
#include <algorithm>

#include "tools.hpp"
#include "problem.hpp"
#include "parameters.hpp"
#include "engine.hpp"
#include "direct_common.hpp"
#include "direct_dict.hpp"
#include "direct_producer.hpp"

/*
 * The direct scheme: the exhaustive meet-in-the-middle over the distributed dictionary (PROTOCOL.md §8).
 * ceil(2^n / w') rounds of two phases each, w' = fill * w entries: FILL inserts f on the round's chunk of
 * the domain, PROBE probes g on the whole domain.  Every evaluation is one point through the engine, so
 * the node runs at its comm thread's speed (PROBLEM.md §2); this is the baseline PCS is measured against.
 * The problem wrappers, the Scheme's functions and the two entry points, claw_search / collision_search.
 */

namespace mitm::direct {

/*
 * A claw problem, f, g : {0,1}^n -> {0,1}^m.  FILL evaluates f, PROBE evaluates g; a match (x, y) is
 * verified with f(x) and tested with is_good_pair(x, y).
 */
template <class Problem>
class ClawWrapper {
public:
	const Problem &pb;
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	const u64 out_mask;             /* the low m bits */
	static constexpr int vlen = Problem::vlen;

	ClawWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
	{
		static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
		              "problem not derived from mitm::AbstractClawProblem");
		assert(n <= 63 && m <= 64);
	}

	/* the FILL phase's function, scalar: what a match is verified with */
	u64 fill(u64 x) const
	{
		return pb.f(x);
	}

	/* vlen evaluations of the phase's function: f while filling, g while probing */
	void veval(u64 phase, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1) {
			y[0] = (phase == FILL) ? pb.f(x[0]) : pb.g(x[0]);
		} else {
			bool choice[vlen];
			for (int k = 0; k < vlen; k++)
				choice[k] = (phase == FILL);
			pb.vfg(x, choice, y);
		}
	}

	/* is the collision f(x) == g(y) the one we want?  The pair comes out in the order the problem wants */
	bool good(u64 &x, u64 &y) const
	{
		return pb.is_good_pair(x, y);
	}

	/* one value every rank must agree on (PROTOCOL.md §4.1): the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.g(b & in_mask);
	}
};

/*
 * A collision problem, f : {0,1}^n -> {0,1}^m: both phases evaluate f.  A match is a pair of distinct
 * preimages; is_good_pair may be ordered, so both orders are tried and the pair is returned in the
 * order that passed.  Unlike PCS's, this collision search works: it is exhaustive.
 */
template <class Problem>
class CollisionWrapper {
public:
	const Problem &pb;
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	const u64 out_mask;             /* the low m bits */
	static constexpr int vlen = Problem::vlen;

	CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
	{
		static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
		              "problem not derived from mitm::AbstractCollisionProblem");
		assert(n <= 63 && m <= 64);
	}

	u64 fill(u64 x) const
	{
		return pb.f(x);
	}

	void veval(u64, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1)
			y[0] = pb.f(x[0]);
		else
			pb.vf(x, y);
	}

	bool good(u64 &x, u64 &y) const
	{
		if (x == y)
			return false;
		if (pb.is_good_pair(x, y))
			return true;
		if (not pb.is_good_pair(y, x))
			return false;
		std::swap(x, y);
		return true;
	}

	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.f(b & in_mask);
	}
};


/******************************** the scheme's functions ********************************/

/* (round, FILL) -> (round, PROBE) -> (round + 1, FILL) ...; past the last round, `stop` (PROTOCOL.md §8) */
template <class Wrapper>
void Scheme::next_header(const Params &params, const Wrapper &, PRNG &, Header &h, u64 &stop)
{
	if (h.phase == PROBE) {
		h.round += 1;
		h.phase = FILL;
	} else {
		h.phase = PROBE;
	}
	if (h.round >= params.n_rounds)
		stop = 1;
}

inline void Scheme::build_dict(SharedContext<Scheme> &shared, const Params &params, int d)
{
	shared.shards[d] = std::make_unique<DirectDict>(params.w_shard, params.n);
}

template <class Wrapper>
void Scheme::producer_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
                             SharedContext<Scheme> &shared, int index)
{
	direct::producer_thread(ctx, wrapper, params, shared, index);
}

template <class Wrapper>
void Scheme::dict_thread(ThreadContext<Scheme> &ctx, const Wrapper &wrapper, const Params &params,
                         SharedContext<Scheme> &shared, int index)
{
	direct::dict_thread(ctx, wrapper, params, shared, index);
}

/* the dictionary lives from FILL to PROBE and is emptied after that (PROTOCOL.md §8) */
inline void Scheme::after_round(SharedContext<Scheme> &shared, const Params &, int d, const Header &h)
{
	if (h.phase == PROBE)
		shared.shards[d]->flush();
}

/* never on volume: a phase ends when every producer has exhausted its range (PROTOCOL.md §4.5) */
inline bool Scheme::round_complete(const Params &, const u64[])
{
	return false;
}

/* the round of protocol round `nround`, and its phase: they alternate from (0, FILL) */
static inline u64 round_of(u64 nround)
{
	return (nround - 1) / 2;
}

static inline const char *phase_name(u64 nround)
{
	return ((nround - 1) % 2 == 0) ? "FILL" : "PROBE";
}

/* the startup report, printed by run() before the team exists */
inline void Scheme::banner(const Params &params, u64 seed)
{
	printf("Starting MPI+OpenMP direct (exhaustive) meet-in-the-middle with seed=%016" PRIx64 "\n", seed);
	layout_banner(params);
	char hper[8], hwork[8];
	human_format(params.per_round, hper);
	human_format((double) params.n_rounds * (params.per_round + params.domain), hwork);
	printf("Dictionary: linear probing, 8-byte slots = %d-bit preimage | %d check bits | occupancy.  "
	       "Fill %.2f: %s entries per round\n", params.n, params.check_bits, params.fill, hper);
	printf("%" PRIu64 " round(s) of two phases: FILL (f on 2^%.2f preimages) then PROBE (g on all 2^%d).  "
	       "Total work: %s evaluations = 2^%.2f\n", params.n_rounds, std::log2((double) params.per_round),
	       params.n, hwork, std::log2((double) params.n_rounds * (params.per_round + params.domain)));
	if (params.capped)
		printf("NOTICE: --nrounds caps the search at %" PRIu64 " round(s): it is NOT exhaustive\n",
		       params.n_rounds);
	printf("NOTE: every evaluation is a point through the comm thread: the node runs at its router's speed"
	       " (PROBLEM.md §2), not its CPUs'.  One rank per L3 domain helps (PROBLEM.md §5.A)\n");
	fflush(stdout);
}

/* the live one-line display, from the reports so far */
inline void Scheme::display(const Params &params, const u64 reported[], double delta, u64 nround)
{
	u64 total = (phase_name(nround)[0] == 'F')
	          ? std::min(params.per_round, params.domain - round_of(nround) * params.per_round)
	          : params.domain;
	double completion = (double) reported[N_POINTS] / total;
	char hrate[8], hdict[8], hnrate[8];
	human_format((double) reported[N_EVAL] / params.n_producers / delta, hrate);
	human_format((double) (reported[N_INSERT] + reported[N_PROBE]) / params.n_dicts / delta, hdict);
	human_format((double) reported[N_POINTS] * POINT_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);
	printf("\rRound %" PRIu64 "/%" PRIu64 " %s:  %.1fs (%.1f%%, ETA %.1fs).  %s #f/s per producer.  "
	       "%s points/s per dict thread.  node-->%sB/s   ",
	       round_of(nround), params.n_rounds, phase_name(nround), delta, 100. * completion,
	       (completion > 0) ? delta * (1 - completion) / completion : 0., hrate, hdict, hnrate);
	fflush(stdout);
}

/* the round report: `r` is the exact SUM of the counters over every thread of every node, `total` the all-time sums */
inline void Scheme::round_report(const Params &params, const u64 r[], const u64 total[], const RoundStats &,
                                 const RoundStats &, double delta, u64 nround)
{
	char hrate[8], hnrate[8], hdict[8];
	human_format((double) r[N_EVAL] / params.n_producers / delta, hrate);
	human_format((double) r[N_POINTS] * POINT_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);
	human_format((double) (r[N_INSERT] + r[N_PROBE]) / params.n_dicts / delta, hdict);
	printf("\n");
	printf("Round %" PRIu64 " %s.  %.1fs.  2^%.2f evaluations (total 2^%.2f).  %s #f/s per producer.  "
	       "node-->%sB/s.  %s points/s per dict thread\n",
	       round_of(nround), phase_name(nround), delta, std::log2((double) (r[N_EVAL] ? r[N_EVAL] : 1)),
	       std::log2((double) (total[N_EVAL] ? total[N_EVAL] : 1)), hrate, hnrate, hdict);
	if (r[N_INSERT] > 0)
		printf("            %" PRIu64 " inserted (%.2f%% of the points shipped), load %.2f/slot, "
		       "%.2f slots visited per insert\n", r[N_INSERT], 100. * r[N_INSERT] / r[N_POINTS],
		       (double) r[N_INSERT] / params.w, (double) r[N_STEPS] / r[N_INSERT]);
	if (r[N_PROBE] > 0)
		printf("            %" PRIu64 " probed (%.2f%% of the points shipped), %.2f slots visited per probe.  "
		       "%" PRIu64 " matches: %" PRIu64 " false positives, %" PRIu64 " collisions (total 2^%.2f)\n",
		       r[N_PROBE], 100. * r[N_PROBE] / r[N_POINTS], (double) r[N_STEPS] / r[N_PROBE], r[N_MATCH],
		       r[BAD_MATCH], r[N_COLLISIONS], std::log2((double) (total[N_COLLISIONS] ? total[N_COLLISIONS] : 1)));
	if (r[STALL_OUT] | r[STALL_IN])
		printf("            STALLED  %" PRIu64 " comm turns holding points / %" PRIu64 " turns with a delivery"
		       " waiting on a full dict ring\n", r[STALL_OUT], r[STALL_IN]);
	printf("\n");
	fflush(stdout);
}

/* the last line: found, or the domain exhausted (a proof of absence, unless --nrounds cut it short) */
inline void Scheme::done(const Params &params, u64 nround, bool found, double seconds)
{
	u64 rounds_done = nround / 2;       /* two protocol rounds per direct round, the stop header not counted */
	if (found)
		printf("Completed in %.2fs\n", seconds);
	else if (params.capped)
		printf("Gave up after %" PRIu64 " of the %" PRIu64 " round(s) an exhaustive search needs (%.2fs)\n",
		       rounds_done, (params.domain + params.per_round - 1) / params.per_round, seconds);
	else
		printf("No solution: the domain is exhausted after %" PRIu64 " round(s) (%.2fs)\n",
		       rounds_done, seconds);
	fflush(stdout);
}


/******************************** the entry points ********************************/

/* find x0 != x1 with f(x0) == f(x1) and is_good_pair(x0, x1), or prove there is none */
template <class Problem>
optional<pair<u64, u64>> collision_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	if (opts.verbose && rank == 0)
		printf("Starting direct collision search with f : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
		       pb.n, pb.m, Problem::vlen);

	CollisionWrapper<Problem> wrapper(pb);
	auto collision = run<Scheme>(wrapper, nbytes_memory, opts, prng);
	if (not collision)
		return nullopt;

	auto [round, x0, x1] = *collision;
	assert(x0 != x1);
	assert(pb.f(x0) == pb.f(x1));
	assert(pb.is_good_pair(x0, x1));
	return optional(pair(x0, x1));
}

/* find x0, x1 with f(x0) == g(x1) and is_good_pair(x0, x1), or prove there is none */
template <class Problem>
optional<pair<u64, u64>> claw_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	if (opts.verbose && rank == 0)
		printf("Starting direct claw search with f, g : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
		       pb.n, pb.m, Problem::vlen);

	ClawWrapper<Problem> wrapper(pb);
	auto claw = run<Scheme>(wrapper, nbytes_memory, opts, prng);
	if (not claw)
		return nullopt;

	auto [round, x0, x1] = *claw;
	assert(pb.f(x0) == pb.g(x1));
	assert(pb.is_good_pair(x0, x1));
	return optional(pair(x0, x1));
}

}
#endif
