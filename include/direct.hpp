#ifndef MITM_DIRECT
#define MITM_DIRECT

#include <mpi.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <memory>
#include <vector>

#include "tools.hpp"
#include "problem.hpp"
#include "parameters.hpp"
#include "router/router.hpp"
#include "direct_common.hpp"
#include "direct_dict.hpp"
#include "direct_producer.hpp"

/*
 * The direct engine: the exhaustive meet-in-the-middle over a distributed dictionary, on the Router.
 * ceil(2^n / w') rounds of two phases each, w' = fill * w entries: FILL inserts f on the round's chunk of the
 * domain, PROBE probes g on the whole domain.  One Router round per phase, so the Router's own end of round
 * is the barrier the algorithm needs -- every entry of the round is in its shard before the first probe of
 * the round is delivered -- and nothing is dropped, so "no solution" is a proof.
 * The problem wrappers, the three thread rounds, the printing and the two entry points, claw_search /
 * collision_search.
 */

namespace mitm::direct {

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
	void veval(int phase, const u64 x[], u64 y[]) const
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

	/* one value every rank must agree on: the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.g(b & in_mask);
	}
};

/*
 * A collision problem, f : {0,1}^n -> {0,1}^m: both phases evaluate f.  A match is a pair of distinct
 * preimages; is_good_pair may be ordered, so both orders are tried and the pair is returned in the order
 * that passed.
 */
template <class Problem>
class CollisionWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	const u64 out_mask;             /* the low m bits */
	static constexpr int vlen = Problem::vlen;

	CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)),
	                                      out_mask(make_mask(pb.m))
	{
		static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
		              "problem not derived from mitm::AbstractCollisionProblem");
		assert(n <= 63 && m <= 64);
	}

	u64 fill(u64 x) const
	{
		return pb.f(x);
	}

	void veval(int, const u64 x[], u64 y[]) const
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


/******************************** the printing ********************************/

/* the startup report, on rank 0 before the team exists; the Router prints its own layout after it */
static void banner(const Params &params, u64 seed)
{
	printf("Starting the direct (exhaustive) meet-in-the-middle on the Router with seed=%016" PRIx64 "\n", seed);
	printf("MPI: %d node(s) x (1 service + %d dict + %d producer) = %d threads/node\n", params.n_nodes,
	       params.R, params.S, params.n_threads);
	printf("MPI: %d dictionary shards, %d producer threads in total\n", params.n_dicts, params.n_producers);
	char hper[8], hwork[8], hslots[8];
	human_format(params.per_round, hper);
	human_format((double) params.n_rounds * (params.per_round + params.domain), hwork);
	human_format(params.w, hslots);
	printf("Dictionary: linear probing, 8-byte slots = %d-bit preimage | %d check bits | occupancy.  "
	       "Fill %.2f: %s entries per round\n", params.n, params.check_bits, params.fill, hper);
	printf("RAM per node == %.1f MB of dictionary; %s slots in all (2^%.2f), %s per shard\n",
	       (double) params.w_shard * params.R * sizeof(u64) / 1e6, hslots, std::log2((double) params.w), hper);
	printf("%" PRIu64 " round(s) of two phases: FILL (f on 2^%.2f preimages) then PROBE (g on all 2^%d).  "
	       "Total work: %s evaluations = 2^%.2f\n", params.n_rounds, std::log2((double) params.per_round),
	       params.n, hwork, std::log2((double) params.n_rounds * (params.per_round + params.domain)));
	if (params.capped)
		printf("NOTICE: --nrounds caps the search at %" PRIu64 " round(s): it is NOT exhaustive\n",
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
	char hrate[8], hdict[8], hnrate[8];
	human_format((double) eval / params.S / delta, hrate);
	human_format((double) retired / params.R / delta, hdict);
	human_format((double) stats[ROUTER_BYTES_SENT] / delta, hnrate);
	printf("\rRound %" PRIu64 "/%" PRIu64 " %s:  %.1fs (%.1f%%, ETA %.1fs).  %s #f/s per producer.  "
	       "%s points/s per dict thread.  node-->%sB/s   ", round, params.n_rounds,
	       (phase == FILL) ? "FILL" : "PROBE", delta, 100. * completion,
	       (completion > 0) ? delta * (1 - completion) / completion : 0., hrate, hdict, hnrate);
	fflush(stdout);
}

/* the end-of-phase report: `r` is the exact sum over every thread of every node, `total` the all-time sums */
static void round_report(const Params &params, const u64 r[], const u64 total[], double delta, u64 round,
                         int phase)
{
	char hrate[8], hnrate[8], hdict[8];
	human_format((double) r[N_EVAL] / params.n_producers / delta, hrate);
	human_format((double) r[REC_ROUTER + ROUTER_BYTES_SENT] / params.n_nodes / delta, hnrate);
	human_format((double) (r[N_INSERT] + r[N_PROBE]) / params.n_dicts / delta, hdict);
	printf("\n");
	printf("Round %" PRIu64 " %s.  %.1fs.  2^%.2f evaluations (total 2^%.2f).  %s #f/s per producer.  "
	       "node-->%sB/s.  %s points/s per dict thread\n", round, (phase == FILL) ? "FILL" : "PROBE", delta,
	       std::log2((double) (r[N_EVAL] ? r[N_EVAL] : 1)),
	       std::log2((double) (total[N_EVAL] ? total[N_EVAL] : 1)), hrate, hnrate, hdict);
	if (r[N_INSERT] > 0)
		printf("            %" PRIu64 " inserted, load %.2f/slot, %.2f slots visited per insert\n",
		       r[N_INSERT], (double) r[N_INSERT] / params.w, (double) r[N_STEPS] / r[N_INSERT]);
	if (r[N_PROBE] > 0)
		printf("            %" PRIu64 " probed, %.2f slots visited per probe.  %" PRIu64 " matches: %" PRIu64
		       " false positives, %" PRIu64 " collisions (total 2^%.2f)\n", r[N_PROBE],
		       (double) r[N_STEPS] / r[N_PROBE], r[N_MATCH], r[BAD_MATCH], r[N_COLLISIONS],
		       std::log2((double) (total[N_COLLISIONS] ? total[N_COLLISIONS] : 1)));
	if (r[REC_ROUTER + ROUTER_STALL_OUT] | r[REC_ROUTER + ROUTER_STALL_IN])
		printf("            STALLED  %" PRIu64 " blocks held back for a full destination / %" PRIu64
		       " receives left unposted for want of a free block\n", r[REC_ROUTER + ROUTER_STALL_OUT],
		       r[REC_ROUTER + ROUTER_STALL_IN]);
	if (r[REC_ROUTER + ROUTER_PUSHED] != r[REC_ROUTER + ROUTER_POPPED])
		printf("            ACCOUNTING BROKEN: %" PRIu64 " points pushed, %" PRIu64 " delivered\n",
		       r[REC_ROUTER + ROUTER_PUSHED], r[REC_ROUTER + ROUTER_POPPED]);
	printf("\n");
	fflush(stdout);
}

/* the last line: found, or the domain exhausted -- a proof of absence, unless --nrounds cut it short */
static void done(const Params &params, u64 phases, bool found, double seconds)
{
	u64 rounds = (phases + 1) / 2;
	if (found)
		printf("Completed in %.2fs\n", seconds);
	else if (params.capped)
		printf("Gave up after %" PRIu64 " of the %" PRIu64 " round(s) an exhaustive search needs (%.2fs)\n",
		       rounds, (params.domain + params.per_round - 1) / params.per_round, seconds);
	else
		printf("No solution: the domain is exhausted after %" PRIu64 " round(s) (%.2fs)\n", rounds, seconds);
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
		rec[REC_X] = shared.golden[1];
		rec[REC_Y] = shared.golden[2];
	}
	MPI_Allgather(rec, REC_WORDS, MPI_UINT64_T, records, REC_WORDS, MPI_UINT64_T, params.mpi_comm);

	u64 sum[REC_WORDS] = {};
	for (int r = 0; r < params.n_nodes; r++) {
		const u64 *q = records + (size_t) r * REC_WORDS;
		for (int k = 0; k < N_COUNTERS + ROUTER_STATS_SIZE; k++)
			sum[k] += q[k];
		if (q[REC_FOUND] && not shared.solved) {
			shared.solved = true;
			shared.solution[0] = round;
			shared.solution[1] = q[REC_X];
			shared.solution[2] = q[REC_Y];
		}
	}
	for (int k = 0; k < N_COUNTERS + ROUTER_STATS_SIZE; k++)
		total[k] += sum[k];
	shared.phases += 1;
	shared.stop = shared.solved || (phase == PROBE && round + 1 >= params.n_rounds);
	if (params.verbose)
		round_report(params, sum, total, delta, round, phase);
}

/* one value every rank must agree on: a mismatch means the ranks hold different problem instances */
template <class Wrapper>
static void self_test(const Wrapper &wrapper, const Params &params, PRNG &prng)
{
	u64 mine[3];
	mine[0] = prng.rand();
	mine[1] = prng.rand();
	mine[2] = wrapper.self_test(mine[0], mine[1]);
	u64 theirs[3] = {mine[0], mine[1], mine[2]};
	MPI_Bcast(theirs, 3, MPI_UINT64_T, 0, params.mpi_comm);
	if (theirs[0] != mine[0] || theirs[1] != mine[1] || theirs[2] != mine[2])
		errx(1, "direct: the ranks do not hold the same problem instance (self-test mismatch)");
}


/******************************** the engine ********************************/

/*
 * The search itself: one OpenMP team per node -- thread 0 the Router's service thread, the next R its dict
 * threads, the rest its producers -- running the same deterministic sequence of phases, (0, FILL), (0, PROBE),
 * (1, FILL), ..., each one Router round.  Returns the round, x and y of the golden pair, the same on every
 * node, or nothing once the domain is exhausted.
 */
template <class Wrapper>
optional<tuple<u64, u64, u64>> run(const Wrapper &wrapper, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "direct: this MPI does not provide MPI_THREAD_FUNNELED");

	const Params params(opts, nbytes_memory, wrapper.n, wrapper.m);
	self_test(wrapper, params, prng);
	if (params.verbose)
		banner(params, prng.seed);

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
		Router_thread rt = Router_Init(role, ROUTER_GROUP_AUTO, params.mpi_comm, ROUTER_TAG, false,
		                               &params.router);
		u64 *ctr = shared.tally[tid].ctr;

		std::unique_ptr<DirectDict> dict;      /* a dict thread's shard, zero-filled here: its own first touch */
		if (role == ROUTER_RECEIVER)
			dict.reset(new DirectDict(params.w_shard, params.n));

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
				dict_round(rt, wrapper, shared, *dict, ctr, round, phase);
			else
				producer_round(rt, wrapper, params, ctr, round, phase);

			Router_Reset(rt);      /* its own team barriers and MPI_Barrier: every tally is written by now */

			if (role == ROUTER_SERVICE)
				epilogue(params, shared, stats.data(), records.data(), total.data(), round, phase,
				         wtime() - t0);

			#pragma omp barrier    /* the verdict, and every thread's next tally, are thread 0's to publish */
			if (shared.stop)
				break;
		}
	}

	if (params.verbose)
		done(params, shared.phases, shared.solved, wtime() - t_start);
	if (not shared.solved)
		return nullopt;
	return optional(tuple(shared.solution[0], shared.solution[1], shared.solution[2]));
}


/******************************** the entry points ********************************/

/* find x0 != x1 with f(x0) == f(x1) and is_good_pair(x0, x1), or prove there is none */
template <class Problem>
optional<pair<u64, u64>> collision_search(const Problem &pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	if (opts.verbose && rank == 0)
		printf("Starting direct collision search with f : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
		       pb.n, pb.m, Problem::vlen);

	CollisionWrapper<Problem> wrapper(pb);
	auto collision = run(wrapper, nbytes_memory, opts, prng);
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
optional<pair<u64, u64>> claw_search(const Problem &pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	if (opts.verbose && rank == 0)
		printf("Starting direct claw search with f, g : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
		       pb.n, pb.m, Problem::vlen);

	ClawWrapper<Problem> wrapper(pb);
	auto claw = run(wrapper, nbytes_memory, opts, prng);
	if (not claw)
		return nullopt;

	auto [round, x0, x1] = *claw;
	assert(pb.f(x0) == pb.g(x1));
	assert(pb.is_good_pair(x0, x1));
	return optional(pair(x0, x1));
}

}
#endif
