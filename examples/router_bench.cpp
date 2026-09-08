#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <vector>
#include <cstdio>
#include <inttypes.h>

#include "tools.hpp"
#include "router/router.hpp"
#include "router_driver.hpp"

/*
 * The Router's throughput benchmark: senders push as fast as they can, one PRNG call per point (a producer
 * computes, it does not idle, and an idling core is clocked down), to uniform, skewed or node-local
 * destinations; receivers pop the points one at a time and hash each (murmur128), folding the hashes by XOR; the
 * service thread
 * runs until quiescent.  Prints, per round, the aggregate rate of points routed (delivered to a receiver), then
 * per node the push and pop rates, the network traffic, the drops, the service thread's duty cycle, the cost of
 * a point on a sender (the PRNG call and the destination pick included) and the XOR of every hash folded over
 * the receivers and the nodes: it keeps the hashes from being optimised away, and in lossless mode it depends
 * on the points alone, not on their order or their route.  Before the rounds, the raw benchmark: the senders
 * run that same loop for a second with Router_Push taken out of it, so the difference between the raw rate and
 * the routed rate is the router and nothing else.
 */
using namespace mitm;

static const u64 seed = 1337;       /* the senders' PRNG streams: this key, Router_rank as the sequence */
std::vector<double> ns_per_point;   /* per sender: a point's cost, the PRNG call and the destination pick included */
u64 xor_hash = 0;                   /* the node's receivers' hashes, folded; the raw run's outputs */
double raw_rate = 0;                /* points/s the node's senders generate without pushing, summed */

/* Where a point goes: the sender's fan-out, and the raw run's, so that the two cannot drift apart.  The modulo
 * is deliberately 32-bit: a 64-bit divq is microcoded on Skylake (36-95 cycles against 26 for divl) and costs
 * more than the push it feeds, which throttles the senders below the router's own rate and hides it.  F and per
 * are tiny, so folding x to 32 bits first biases nothing that matters here. */
struct Fanout {
	int F;                      /* destinations: the receivers, or --dests' virtual fan-out */
	int per;                    /* local-only: destinations on this node, [base, base + per) */
	int base;
	bool local_only;

	Fanout(const Router_thread &rt, const RouterArgs &a)
	    : F(rt.node.F), per(rt.node.per_node), base(rt.node.rank * rt.node.per_node), local_only(a.local_only)
	{
	}

	int operator()(u64 x) const
	{
		return local_only ? base + (int) ((u32) x % (u32) per) : (int) ((u32) x % (u32) F);
	}
};

/* a sender runs the sending loop without Router_Push for a second: the generation rate the routing is measured
 * against.  Everything the sender does per point but the push itself, the destination pick included */
static void bench_raw(Router_thread &rt, const RouterArgs &a)
{
	Fanout fan(rt, a);
	PRNG prng(seed, (u64) Router_rank(rt));
	u64 acc = 0;
	u64 n = 0;
	double t0 = wtime();
	double t;
	for (;;) {
		for (int j = 0; j < 1024; j++) {
			u64 x = prng.rand();
			acc ^= x ^ (u64) fan(x);   /* consumed, or the destination pick is optimised away */
		}
		n += 1024;
		t = wtime();
		if (t >= t0 + 1)
			break;
	}
	#pragma omp atomic
	raw_rate += (double) n / (t - t0);
	#pragma omp atomic
	xor_hash ^= acc;
}

/* rank 0 prints the raw benchmark: the senders' aggregate generation rate over all nodes */
static void raw_report(const Router_thread &rt)
{
	double total;
	u64 folded;
	MPI_Reduce(&raw_rate, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&xor_hash, &folded, 1, MPI_UINT64_T, MPI_BXOR, 0, MPI_COMM_WORLD);
	if (rt.node.rank != 0)
		return;
	int nsend = Router_num_send(rt);
	char raw_s[8];
	human_format((u64) total, raw_s);
	printf("raw: %d senders generate %s pts/s without pushing (%.1f ns/point) | xor %016" PRIx64 "\n",
	       nsend, raw_s, 1e9 * nsend / total, folded);
}

static void bench_sender(Router_thread &rt, const RouterArgs &a, double t_end)
{
	Fanout fan(rt, a);
	PRNG prng(seed, (u64) Router_rank(rt));
	u64 n = 0;
	double t0 = wtime();
	for (;;) {
		for (int j = 0; j < 1024; j++) {
			u64 x = prng.rand();
			Router_Push(x, n, fan(x), rt);
			n += 1;
		}
		if (wtime() >= t_end)
			break;
	}
	double busy = wtime() - t0;
	Router_Close(rt);
	ns_per_point[Router_local_rank(rt)] = 1e9 * busy / (double) n;
}

static void bench_receiver(Router_thread &rt)
{
	u64 x = 0;
	for (;;) {
		u64 a, b;
		if (not Router_Pop(&a, &b, rt)) {
			if (Router_Test_drained(rt))
				break;
			cpu_relax();
			continue;
		}
		x ^= murmur128(a, b);
	}
	#pragma omp atomic
	xor_hash ^= x;
}

/* rank 0 prints the round: the aggregate routed rate, sums over the nodes, the per-node spread of the push
 * rate, the folded hash */
static void bench_report(const Router_thread &rt, const RouterArgs &a, int round, double elapsed)
{
	u64 st[ROUTER_STATS_SIZE];
	Router_Stats(st, rt);
	double mine[4] = {(double) st[ROUTER_PUSHED] / elapsed, 0, 0, elapsed};
	for (int s = 0; s < a.senders; s++)
		mine[1] += ns_per_point[s] / a.senders;
	mine[2] = (double) (st[ROUTER_TURNS] - st[ROUTER_IDLE_TURNS]) / (double) (st[ROUTER_TURNS] ? st[ROUTER_TURNS] : 1);
	u64 tot[ROUTER_STATS_SIZE];
	u64 folded;                         /* the nodes' xor_hash, folded */
	double lo[4];
	double hi[4];
	double sum[4];
	MPI_Reduce(st, tot, ROUTER_STATS_SIZE, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&xor_hash, &folded, 1, MPI_UINT64_T, MPI_BXOR, 0, MPI_COMM_WORLD);
	MPI_Reduce(mine, lo, 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(mine, hi, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(mine, sum, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rt.node.rank != 0)
		return;
	int P = rt.node.n_nodes;
	double t = hi[3];
	u64 drops = tot[ROUTER_DROPPED_SERVICE] + tot[ROUTER_DROPPED_NET] + tot[ROUTER_DROPPED_RECV];
	char routed_s[8];
	char push_s[8];
	char pop_s[8];
	char net_s[8];
	char msg_s[8];
	char blk_s[8];
	char turn_s[8];
	human_format((u64) (tot[ROUTER_POPPED] / t), routed_s);
	human_format((u64) (tot[ROUTER_TURNS] / t / P), turn_s);
	human_format((u64) (tot[ROUTER_PUSHED] / t / P), push_s);
	human_format((u64) (tot[ROUTER_POPPED] / t / P), pop_s);
	human_format((u64) (tot[ROUTER_SENT] / t / P), net_s);
	human_format((u64) (tot[ROUTER_MSGS_SENT] / t / P), msg_s);
	human_format((u64) (tot[ROUTER_BLOCKS] / t / P), blk_s);
	printf("round %d: %.2fs | routed %s pts/s | per node: push %s/s (%.0f-%.0f M/s) pop %s/s net %s pts/s %.2f GB/s"
	       " %s msgs/s %s blocks/s | %.1f ns/point | service %.0f%% busy (%s turns/s)"
	       " | dropped %.3f%% (%" PRIu64 " svc, %" PRIu64 " net, %" PRIu64 " recv) | xor %016" PRIx64 "\n",
	       round, t, routed_s, push_s, lo[0] / 1e6, hi[0] / 1e6, pop_s, net_s,
	       (double) tot[ROUTER_BYTES_SENT] / t / P / 1e9, msg_s, blk_s, sum[1] / P,
	       100. * sum[2] / P, turn_s,
	       tot[ROUTER_PUSHED] ? 100. * drops / tot[ROUTER_PUSHED] : 0.,
	       tot[ROUTER_DROPPED_SERVICE], tot[ROUTER_DROPPED_NET], tot[ROUTER_DROPPED_RECV], folded);
	if (tot[ROUTER_PUSHED] != tot[ROUTER_POPPED] + drops)
		printf("  ACCOUNTING BROKEN: pushed %" PRIu64 " != popped %" PRIu64 " + dropped %" PRIu64 "\n",
		       tot[ROUTER_PUSHED], tot[ROUTER_POPPED], drops);
}

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI_THREAD_FUNNELED not provided");
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	RouterArgs a;
	router_parse(argc, argv, a);
	a.opts.verbose = a.opts.verbose && rank == 0;
	int S = a.senders;
	int R = a.receivers;

	ns_per_point.assign(S, 0);
	double t0 = 0;

	#pragma omp parallel num_threads(1 + S + R)
	{
		int tid = omp_get_thread_num();
		int role = (tid == 0) ? ROUTER_SERVICE : (tid <= R ? ROUTER_RECEIVER : ROUTER_SENDER);
		Router_thread rt = Router_Init(role, ROUTER_GROUP_AUTO, MPI_COMM_WORLD, 42, a.lossy, &a.opts);
		if (role == ROUTER_SENDER)
			bench_raw(rt, a);
		#pragma omp barrier
		if (tid == 0)
			raw_report(rt);
		for (int round = 0; round < a.rounds; round++) {
			if (tid == 0) {
				xor_hash = 0;
				t0 = wtime();
			}
			#pragma omp barrier
			if (tid == 0) {
				while (not Router_Test_quiescent(rt))
					Router_Progress(rt);
			} else if (role == ROUTER_RECEIVER) {
				bench_receiver(rt);
			} else {
				bench_sender(rt, a, t0 + a.seconds);
			}
			#pragma omp barrier
			if (tid == 0) {
				bench_report(rt, a, round, wtime() - t0);
				Router_Reset(rt);
			}
		}
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}
