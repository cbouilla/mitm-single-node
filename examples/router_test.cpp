#include <mpi.h>
#include <omp.h>
#include <atomic>
#include <vector>
#include <cstdio>
#include <cstring>
#include <inttypes.h>

#include "router/router.hpp"
#include "router_driver.hpp"

/*
 * The Router's test suite.  Every test is one configuration of the same round: senders push
 * a deterministic stream and close, receivers read blocks in place and verify until drained, the service thread runs
 * Router_Progress until quiescent, then rank 0 checks the global accounting.  CHECK, never assert: the
 * build may carry -DNDEBUG.
 */
using namespace mitm;

static int g_rank = 0;
static int g_nodes = 1;
static std::atomic<int> g_failed(0);

#define CHECK(cond) do { if (!(cond)) { \
	fprintf(stderr, "CHECK failed on rank %d: %s (%s:%d)\n", g_rank, #cond, __FILE__, __LINE__); \
	g_failed.fetch_add(1); } } while (0)

static const u64 SALT = 0x5a17c0ffee5a17ull;

/* the stream of sender s of node `node` in round r: a encodes (r, node, s, i), b encodes a and the dest */
static inline u64 point_a(int round, int node, int s, u64 i)
{
	return ((u64) round << 56) | ((u64) node << 48) | ((u64) s << 40) | i;
}

static inline u64 point_b(u64 a, int dest)
{
	return murmur64(a) ^ SALT ^ (u64) dest;
}

struct Shared {
	std::vector<u64> sent_to;        /* S * F */
	std::vector<u64> got;            /* per local receiver */
	std::vector<u64> seen;           /* per local receiver: a bitmap over (node, s, i) */
	u64 max_points = 0;
	std::vector<int> t_group;        /* per thread: its Router_group, gathered after Init */
	std::vector<int> t_cpu;          /* per thread: its Router_cpu */
	std::vector<int> t_domain;       /* per thread: its Router_domain */
	std::vector<int> t_grcv;         /* per sender thread: its group's first receiver, or -1 */
};

static u64 points_of(const RouterArgs &a, int s)
{
	if (a.stagger)
		return (a.points * (u64) (s + 1)) / (u64) a.senders;
	return a.points;
}

static void sender_round(Router_thread &rt, const RouterArgs &a, Shared &sh, int s, int round)
{
	int F = Router_num_recv(rt);
	u64 *sent_to = sh.sent_to.data() + (size_t) s * F;
	u64 i = 0;
	if (a.partial) {
		for (int d = 0; d < F; d++) {
			int n = 1 + (d + s + round) % 7;
			for (int j = 0; j < n; j++, i++) {
				u64 x = point_a(round, g_rank, s, i);
				sent_to[d] += 1;
				Router_Push(x, point_b(x, d), d, rt);
			}
		}
	} else {
		u64 n = points_of(a, s);
		for (; i < n; i++) {
			u64 x = point_a(round, g_rank, s, i);
			int d = a.skew > 0 ? (int) (murmur64(x) % (u64) a.skew) : (int) (murmur64(x) % (u64) F);
			sent_to[d] += 1;
			Router_Push(x, point_b(x, d), d, rt);
		}
	}
	Router_Close(rt);
}

/* one delivered point: from this round, from a real sender, to this receiver's rank, and never seen before */
static void check_point(u64 x, u64 y, int me, int S, Shared &sh, u64 *seen, int round)
{
	int rnd = (int) (x >> 56);
	int node = (int) ((x >> 48) & 0xff);
	int s = (int) ((x >> 40) & 0xff);
	u64 i = x & ((1ull << 40) - 1);
	CHECK(rnd == round);
	CHECK(y == point_b(x, me));
	CHECK(node < g_nodes && s < S && i < sh.max_points);
	u64 bit = ((u64) (node * S + s)) * sh.max_points + i;
	CHECK((seen[bit / 64] & (1ull << (bit % 64))) == 0);
	seen[bit / 64] |= 1ull << (bit % 64);
}

static void receiver_round(Router_thread &rt, const RouterArgs &a, Shared &sh, int r, int round)
{
	int me = Router_rank(rt);
	int S = a.senders;
	u64 *seen = sh.seen.data() + (size_t) r * ((g_nodes * S * sh.max_points + 63) / 64);
	u64 got = 0;
	for (;;) {
		size_t k = 0;
		if (a.pop_one) {
			u64 x, y;
			if (Router_Pop(&x, &y, rt)) {
				check_point(x, y, me, S, sh, seen, round);
				k = 1;
			}
		} else {
			const u64 *pts;
			k = Router_Grab(&pts, rt);
			for (size_t j = 0; j < k; j++)
				check_point(pts[2 * j], pts[2 * j + 1], me, S, sh, seen, round);
			if (k > 0)
				Router_Release(rt);       /* after the checks: they read block memory */
		}
		if (k == 0) {
			if (Router_Test_drained(rt))
				break;
			cpu_relax();
			continue;
		}
		got += k;
		if (a.slow_recv)
			for (int t = 0; t < 2000; t++)
				cpu_relax();
	}
	sh.got[r] = got;
}

/* rank 0's verdict on the round: every point pushed was popped, and what was sent was received */
static void check_round(const Router_thread &rt, const RouterArgs &a, Shared &sh, int round)
{
	int F = Router_num_recv(rt);
	int S = a.senders;
	int R = a.receivers;
	u64 st[ROUTER_STATS_SIZE];
	Router_Stats(st, rt);

	std::vector<u64> node_sent(F, 0);
	for (int s = 0; s < S; s++)
		for (int d = 0; d < F; d++)
			node_sent[d] += sh.sent_to[(size_t) s * F + d];
	std::vector<u64> node_got(F, 0);
	u64 local_pushed = 0;
	u64 local_popped = 0;
	for (int d = 0; d < F; d++)
		local_pushed += node_sent[d];
	for (int r = 0; r < R; r++) {
		node_got[g_rank * R + r] = sh.got[r];
		local_popped += sh.got[r];
	}
	CHECK(st[ROUTER_PUSHED] == local_pushed);
	CHECK(st[ROUTER_POPPED] == local_popped);

	std::vector<u64> tot_sent(F, 0);
	std::vector<u64> tot_got(F, 0);
	u64 tot[ROUTER_STATS_SIZE];
	MPI_Reduce(node_sent.data(), tot_sent.data(), F, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(node_got.data(), tot_got.data(), F, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(st, tot, ROUTER_STATS_SIZE, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	if (g_rank != 0)
		return;
	CHECK(tot[ROUTER_PUSHED] == tot[ROUTER_POPPED]);
	CHECK(tot[ROUTER_SENT] == tot[ROUTER_RECV]);
	CHECK(tot[ROUTER_MSGS_SENT] == tot[ROUTER_MSGS_RECV]);
	CHECK(tot[ROUTER_BYTES_SENT] == tot[ROUTER_BYTES_RECV]);
	for (int d = 0; d < F; d++)
		CHECK(tot_sent[d] == tot_got[d]);
	if (a.opts.verbose)
		printf("  round %d: pushed %" PRIu64 " popped %" PRIu64 " local %" PRIu64 " net %" PRIu64
		       " msgs %" PRIu64 " blocks %" PRIu64 " turns %" PRIu64 " (%" PRIu64 " idle)\n",
		       round, tot[ROUTER_PUSHED], tot[ROUTER_POPPED], tot[ROUTER_LOCAL], tot[ROUTER_SENT],
		       tot[ROUTER_MSGS_SENT], tot[ROUTER_BLOCKS], tot[ROUTER_TURNS], tot[ROUTER_IDLE_TURNS]);
}

/* rank-local, on thread 0 after Init: the plan the Router built for this node -- the group ids, the pinning,
 * and that each producer's group holds the receiver it will pair with. */
static void check_placement(const Router_thread &rt, const RouterArgs &a, Shared &sh, bool use_colors)
{
	int S = a.senders;
	int R = a.receivers;
	int nt = 1 + S + R;
	int G = Router_num_groups(rt);
	int D = Router_num_domains(rt);
	CHECK(G >= 1);
	CHECK(sh.t_group[0] == -1);                 /* the service has no group */
	std::vector<int> gs(G, 0);
	std::vector<int> gr(G, 0);
	for (int t = 1; t < nt; t++) {
		CHECK(sh.t_group[t] >= 0 && sh.t_group[t] < G);
		if (sh.t_group[t] < 0 || sh.t_group[t] >= G)
			continue;
		if (t <= R)
			gr[sh.t_group[t]] += 1;
		else
			gs[sh.t_group[t]] += 1;
	}
	if (a.opts.pin) {
		for (int t = 0; t < nt; t++) {
			CHECK(sh.t_domain[t] >= 0 && sh.t_domain[t] < D);
			for (int u = t + 1; u < nt; u++)
				CHECK(sh.t_cpu[t] != sh.t_cpu[u]);   /* every thread on a CPU of its own */
		}
	} else {
		CHECK(D == 0);
		for (int t = 0; t < nt; t++)
			CHECK(sh.t_domain[t] == -1);
	}
	if (use_colors) {
		CHECK(G == ((S >= 2 && R >= 2) ? 2 : 1));
	} else if (S >= G && R >= G) {
		for (int g = 0; g < G; g++) {               /* enough of each role for one of each per group */
			CHECK(gs[g] >= 1);
			CHECK(gr[g] >= 1);
		}
		for (int t = R + 1; t < nt; t++) {          /* a sender's group's first receiver is in its group */
			int rcv = sh.t_grcv[t];
			if (rcv >= 0)
				CHECK(sh.t_group[1 + rcv] == sh.t_group[t]);
		}
	}
}

/* one configuration, start to finish: a team, its Router, `rounds` rounds; everything dies with the region */
static void run_config(const RouterArgs &a, const char *name, bool use_colors = false)
{
	if (g_rank == 0 && a.opts.verbose)
		printf("== %s (%d senders, %d receivers, %" PRIu64 " points, %d rounds)\n", name,
		       a.senders, a.receivers, a.points, a.rounds);
	int S = a.senders;
	int R = a.receivers;
	int nt = 1 + S + R;
	Shared sh;
	sh.max_points = a.partial ? (u64) 7 * R * g_nodes : (a.points > 0 ? a.points : 1);
	sh.sent_to.assign((size_t) S * R * g_nodes, 0);
	sh.got.assign(R, 0);
	sh.seen.assign((size_t) R * ((g_nodes * S * sh.max_points + 63) / 64), 0);
	sh.t_group.assign(nt, 0);
	sh.t_cpu.assign(nt, 0);
	sh.t_domain.assign(nt, 0);
	sh.t_grcv.assign(nt, -1);

	#pragma omp parallel num_threads(nt)
	{
		int tid = omp_get_thread_num();
		int role = (tid == 0) ? ROUTER_SERVICE : (tid <= R ? ROUTER_RECEIVER : ROUTER_SENDER);
		int li = (role == ROUTER_RECEIVER) ? (tid - 1) : (tid - 1 - R);   /* the caller's color, in colors mode */
		int group = ROUTER_GROUP_AUTO;
		if (use_colors && role != ROUTER_SERVICE)
			group = (S >= 2 && R >= 2) ? (li % 2) : 0;
		Router_thread rt = Router_Init(role, group, MPI_COMM_WORLD, 42, &a.opts);
		sh.t_group[tid] = Router_group(rt);
		sh.t_cpu[tid] = Router_cpu(rt);
		sh.t_domain[tid] = Router_domain(rt);
		if (role == ROUTER_SENDER && Router_group_num_receivers(rt) > 0)
			sh.t_grcv[tid] = Router_group_receiver(0, rt);
		#pragma omp barrier
		if (tid == 0)
			check_placement(rt, a, sh, use_colors);
		#pragma omp barrier
		for (int round = 0; round < a.rounds; round++) {
			if (tid == 0) {
				for (size_t k = 0; k < sh.sent_to.size(); k++)
					sh.sent_to[k] = 0;
				for (size_t k = 0; k < sh.seen.size(); k++)
					sh.seen[k] = 0;
			}
			#pragma omp barrier
			if (tid == 0) {
				double t0 = wtime();
				while (not Router_Test_quiescent(rt)) {
					Router_Progress(rt);
					if (wtime() - t0 > 10) {
						Router_Dump(stderr, rt);
						MPI_Abort(MPI_COMM_WORLD, 2);
					}
				}
			} else if (role == ROUTER_RECEIVER) {
				receiver_round(rt, a, sh, tid - 1, round);
			} else {
				sender_round(rt, a, sh, tid - 1 - R, round);
			}
			#pragma omp barrier               /* the test's own: check_round reads what the workers wrote to sh */
			if (tid == 0)
				check_round(rt, a, sh, round);
			Router_Reset(rt);
		}
	}
}

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI_THREAD_FUNNELED not provided");
	MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &g_nodes);
	RouterArgs a;
	router_parse(argc, argv, a);
	if (a.opts.dests_per_node != 0)
		errx(1, "router_test: --dests is for the bench, the checks assume real destinations");
	a.opts.verbose = a.opts.verbose && g_rank == 0;
	bool all = a.test == "all";

	if (all || a.test == "connect_only") {
		RouterArgs c = a;
		c.rounds = 0;
		run_config(c, "connect_only");
	}
	if (all || a.test == "basic") {
		RouterArgs c = a;
		run_config(c, "basic");
	}
	if (all || a.test == "tiny") {
		RouterArgs c = a;
		c.opts.block_points = 8;
		c.opts.swc_linesize = 4;
		c.opts.n_recv = 2;
		c.opts.inbox_blocks = 4;
		c.opts.sweep_blocks = 2;
		c.points = a.points / 5;
		run_config(c, "tiny");
	}
	if (all || a.test == "rounds") {
		RouterArgs c = a;
		c.rounds = 5;
		run_config(c, "rounds");
	}
	if (all || a.test == "zero") {
		RouterArgs c = a;
		c.points = 0;
		run_config(c, "zero");
	}
	if (all || a.test == "skew") {
		RouterArgs c = a;
		c.skew = 1;
		run_config(c, "skew");
	}
	if (all || a.test == "partial_lines") {
		RouterArgs c = a;
		c.partial = true;
		c.opts.swc_linesize = 32;
		run_config(c, "partial_lines");
	}
	if (all || a.test == "pop_one") {
		RouterArgs c = a;
		c.pop_one = true;
		run_config(c, "pop_one");
	}
	if (all || a.test == "swc_is_block") {
		RouterArgs c = a;
		c.opts.swc_linesize = 64;
		c.opts.block_points = 64;
		run_config(c, "swc_is_block");
	}
	if (all || a.test == "staggered_close") {
		RouterArgs c = a;
		c.stagger = true;
		run_config(c, "staggered_close");
	}
	if (all || a.test == "slow_recv") {
		RouterArgs c = a;
		c.slow_recv = true;
		c.points = a.points / 20;
		run_config(c, "slow_recv");
	}
	if (all || a.test == "groups") {
		RouterArgs c = a;
		c.opts.pin = false;             /* one L3 domain on a laptop: exercise several groups without pinning */
		c.opts.group_size = 2;
		run_config(c, "groups");
	}
	if (all || a.test == "colors") {
		RouterArgs c = a;
		run_config(c, "colors", true);
	}
	if ((all && g_nodes == 1) || a.test == "asymmetric") {
		RouterArgs c = a;
		c.senders = 1;
		c.receivers = 6;
		run_config(c, "asymmetric");
		c.senders = 6;
		c.receivers = 1;
		run_config(c, "asymmetric");
	}

	int failed = g_failed.load();
	int total = 0;
	MPI_Reduce(&failed, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (g_rank == 0)
		printf(total == 0 ? "router_test: all checks passed\n" : "router_test: %d CHECKS FAILED\n", total);
	MPI_Finalize();
	return total == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
