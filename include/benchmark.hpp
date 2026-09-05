#ifndef MITM_BENCHMARK
#define MITM_BENCHMARK

#include <mpi.h>
#include <cmath>
#include <atomic>
#include <memory>
#include <vector>

#include <omp.h>

#include "tools.hpp"
#include "parameters.hpp"
#include "inserter.hpp"

/*
 * How fast do a problem's f / g (and vfg) iterate, per rank and over the ranks?  For the *_bench
 * drivers; needs no RAM budget, only the Options.
 */

namespace mitm {

/* the rate of N * vlen evaluations since `start`: min, max, mean and std over the ranks, printed by
   rank 0.  Collective */
static void display_stats(u64 N, double start, int vlen, MPI_Comm comm, int rank, int n_nodes)
{
	double rate = vlen * N / (wtime() - start);
	double rate_min = rate;
	double rate_max = rate;
	double rate_avg = rate;
	MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, comm);
	MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, comm);
	rate_avg /= n_nodes;
	double rate_std = (rate - rate_avg) * (rate - rate_avg);
	MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, comm);
	rate_std /= n_nodes;
	rate_std = std::sqrt(rate_std);
	if (rank == 0) {
		char hmin[8], hmax[8], havg[8], hstd[8];
		human_format(rate_min, hmin);
		human_format(rate_max, hmax);
		human_format(rate_avg, havg);
		human_format(rate_std, hstd);
		printf("Benchmark. f/s: min %s max %s avg %s std %s\n", hmin, hmax, havg, hstd);
	}
}

/*
 * Iterate f / g 2^26 times, then vfg 2^20 times if vlen > 1, on every rank, and print the rates.
 * Both loops walk ONE dependent chain per lane, the shape a trail and a collision resolution both
 * have, so that the two rates can be divided: what that ratio is worth is how much a producer gains by
 * batching its resolutions (PROTOCOL.md §3.2).  Each loop's last value is printed, and that is the
 * only reason it survives -O3: nothing else here is observable.  Collective
 */
template<typename Problem>
void benchmark(const Problem& pb, const Options &opts)
{
	int rank, n_nodes;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	MPI_Comm_size(opts.mpi_comm, &n_nodes);

	if (rank == 0)
		printf("Benchmarking scalar implementation (using %d processes)\n", n_nodes);

	MPI_Barrier(opts.mpi_comm);

	u64 mask = make_mask(pb.n);
	u64 N = 1ull << 26;
	u64 x1 = 1;
	double start = wtime();
	for (u64 i = 0; i < N; i++)
		x1 = ((i & 1) ? pb.f(x1) : pb.g(x1)) & mask;
	display_stats(N, start, 1, opts.mpi_comm, rank, n_nodes);
	if (rank == 0)
		printf("  checksum %016" PRIx64 "\n", x1);

	constexpr int vlen = Problem::vlen;
	if constexpr (vlen > 1) {
		if (rank == 0)
			printf("Benchmarking vector implementation (vlen=%d)\n", vlen);

		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 z[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		bool choice[vlen];
		for (int i = 0; i < vlen; i++) {
			choice[i] = i & 1;
			x[i] = i;
		}

		MPI_Barrier(opts.mpi_comm);

		double start = wtime();
		u64 N = 1ull << 20;
		for (u64 i = 0; i < N; i++) {
			pb.vfg(x, choice, z);
			for (int j = 0; j < vlen; j++)
				x[j] = z[j] & mask;
		}
		display_stats(N, start, vlen, opts.mpi_comm, rank, n_nodes);
		u64 check = 0;
		for (int j = 0; j < vlen; j++)
			check ^= x[j];
		if (rank == 0)
			printf("  checksum %016" PRIx64 "\n", check);
	}
}


/*
 * How fast can a dict thread probe its shard, with nothing else running?  R_p has only ever been
 * measured through the whole pipeline, tangled up with the queues and the router (PROBLEM.md §6).
 * This is the shard alone, at the size and the placement the engine would give it: the rank's
 * `dicts_per_node` threads, pinned exactly where Parameters puts them, each hammering its own
 * PcsDict with pseudo-random endpoints.  The rate is reported per quarter-`w` batch, because the
 * probe's cost depends on how full the shard is: a round starts empty and ends at `beta`
 * insertions per slot.  Rank-local, no MPI, no producers, no comm thread.
 */
template<typename Problem>
void probe_benchmark(const Problem& pb, const Options &opts, u64 nbytes_memory)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	Parameters params(opts, nbytes_memory, pb.n, pb.m);
	if (rank != 0)
		return;                     /* every rank is doing the same thing; one of them reports */

	int I = params.dicts_per_node;
	char hw[8], hws[8];
	human_format(params.w, hw);
	human_format(params.w_shard, hws);
	printf("Benchmarking dictionary probes: %d shard(s) of %s slots (%s total, %.1f MB each)\n",
	       I, hws, hw, params.w_shard * sizeof(u64) / 1e6);

	u64 batch = params.w_shard / 4;                   /* probes per timed batch, per shard */
	int n_batches = std::max(1, (int) (4 * params.beta));   /* --beta bounds the run: it is long */
	u64 key_range = std::ldexp(1., pb.n) / params.n_dicts;
	std::vector<double> rate(n_batches * I);
	std::vector<u64> hits(I, 0);
	cpu_set_t caller;                                 /* thread 0 of the team is the calling thread:
	                                                     give it back its mask, or the next benchmark
	                                                     derives its Parameters from one CPU */
	CPU_ZERO(&caller);
	if (sched_getaffinity(0, sizeof(caller), &caller) != 0)
		err(1, "probe_benchmark: sched_getaffinity");

#pragma omp parallel num_threads(I)
	{
		int tid = omp_get_thread_num();
		if (params.bind_threads && pin_to_cpu(params.place.thread_cpu[1 + tid]) < 0)
			warn("probe_benchmark: cannot pin thread %d to CPU %d", tid, params.place.thread_cpu[1 + tid]);
		PcsDict dict(params.jbits, params.w_shard);   /* zero-filled here: NUMA first touch */
		u64 nhit = 0;
		u64 i = 0;
		for (int b = 0; b < n_batches; b++) {
#pragma omp barrier
			double start = wtime();
			for (u64 k = 0; k < batch; k++, i++) {
				u64 end = murmur64(i * 0x9e3779b97f4a7c15ull + tid) % key_range;
				if (dict.pop_insert(end, i, 42))
					nhit += 1;
			}
			rate[b * I + tid] = batch / (wtime() - start);
		}
		hits[tid] = nhit;
	}

	for (int b = 0; b < n_batches; b++) {
		double sum = 0, min = rate[b * I], max = rate[b * I];
		for (int t = 0; t < I; t++) {
			sum += rate[b * I + t];
			min = std::min(min, rate[b * I + t]);
			max = std::max(max, rate[b * I + t]);
		}
		char havg[8], hmin[8], hmax[8], hsum[8];
		human_format(sum / I, havg);
		human_format(min, hmin);
		human_format(max, hmax);
		human_format(sum, hsum);
		printf("  load %.2f -> %.2f/slot: %s probe/s per shard (min %s, max %s), %s over the node\n",
		       0.25 * b, 0.25 * (b + 1), havg, hmin, hmax, hsum);
	}
	if (sched_setaffinity(0, sizeof(caller), &caller) != 0)
		err(1, "probe_benchmark: sched_setaffinity");

	u64 total_hits = 0;
	for (int t = 0; t < I; t++)
		total_hits += hits[t];
	printf("  %" PRIu64 " hits in %" PRIu64 " probes (%.2f%%)\n", total_hits,
	       (u64) n_batches * batch * I, 100. * total_hits / ((double) n_batches * batch * I));
}


/******************************** staging benchmark ***************************/

/* One destination's cursor into its shared staging buffer, alone on its cache line. */
struct alignas(64) StagingCursor {
	std::atomic<u64> next;              /* points ever reserved; the buffer is circular, so it never fills */
};

static const u64 STAGING_POINTS = 1 << 21;   /* points staged per thread and per measurement */
static const u64 STAGING_LINE = 8;           /* u64 in one destination's write-combining line */

/*
 * The floor: a thread appending points to a buffer of its own, with no destination fan-out at all --
 * what the two kernels below pay their fan-out against.  Returns a checksum, the only reason any of
 * these three survives -O3.
 */
template<int WORDS>
static u64 stage_local(u64 *buf, u64 cap, u64 n, u64 seed)
{
	u64 check = 0;
	for (u64 i = 0; i < n; i++) {
		u64 x = murmur64(seed + i);
		u64 *p = buf + (i & (cap - 1)) * WORDS;
		for (int k = 0; k < WORDS; k++)
			p[k] = x + k;
		check += x;
	}
	return check;
}

/*
 * Every producer writing straight into its destination's shared buffer: one atomic reservation and one
 * scattered store per point.  The destination is masked out of the point rather than divided out of
 * it, so that this measures staging and not the runtime division the engine still owes (PROBLEM.md §5.C).
 */
template<int WORDS>
static u64 stage_direct(u64 *buf, StagingCursor *cur, u64 dmask, u64 cap, u64 n, u64 seed)
{
	u64 check = 0;
	for (u64 i = 0; i < n; i++) {
		u64 x = murmur64(seed + i);
		u64 d = x & dmask;
		u64 off = cur[d].next.fetch_add(1, std::memory_order_relaxed) & (cap - 1);
		u64 *p = buf + (d * cap + off) * WORDS;
		for (int k = 0; k < WORDS; k++)
			p[k] = x + k;
		check += off;
	}
	return check;
}

/*
 * The same through one private cache line per destination: the shared cursor and the shared buffer are
 * touched once per full line instead of once per point, at the price of a private line touched on
 * every point.  The fill counts live outside the lines, compact enough to stay in L1.  `cap` and the
 * run are both powers of two and the cursor starts at zero, so a reserved run never straddles the end
 * of a destination's buffer: that is what lets the copy below be one unconditional loop.
 */
template<int WORDS>
static u64 stage_line(u64 *buf, StagingCursor *cur, u64 dmask, u64 cap, u64 n, u64 seed,
                      u64 *wc, u8 *fill)
{
	const u64 run = STAGING_LINE / WORDS;    /* points per line: 2 at three words, 4 at two */
	u64 check = 0;
	for (u64 i = 0; i < n; i++) {
		u64 x = murmur64(seed + i);
		u64 d = x & dmask;
		u64 *l = wc + d * STAGING_LINE;
		u64 f = fill[d];
		for (int k = 0; k < WORDS; k++)
			l[f * WORDS + k] = x + k;
		f += 1;
		if (f < run) {
			fill[d] = f;
			continue;
		}
		fill[d] = 0;
		u64 off = cur[d].next.fetch_add(run, std::memory_order_relaxed) & (cap - 1);
		u64 *p = buf + (d * cap + off) * WORDS;
		for (u64 k = 0; k < run * WORDS; k++)
			p[k] = l[k];
		check += off;
	}
	return check;
}

/*
 * One row of the sweep: `n_dest` destinations, the three strategies in turn, on the rank's producer
 * threads pinned where Parameters puts them.  The shared buffers are first-touched by destination, so
 * they spread over the NUMA nodes the producers sit on -- no producer has a local set of destinations,
 * which is the truth of an all-to-all keyed on a hash of the endpoint.
 */
template<int WORDS>
static void staging_row(const Parameters &params, u64 n_dest, u64 cap, u64 &checksum)
{
	int W = params.producers_per_node;
	u64 dmask = n_dest - 1;
	std::unique_ptr<u64[]> buf(new u64[n_dest * cap * WORDS]);      /* uninitialized: the threads touch it */
	std::unique_ptr<StagingCursor[]> cur(new StagingCursor[n_dest]);
	std::vector<u64> check(W, 0);
	double rate[3] = {};
	double t0 = 0;
	cpu_set_t caller;                   /* thread 0 of the team is the calling thread: give it its mask back */
	CPU_ZERO(&caller);
	if (sched_getaffinity(0, sizeof(caller), &caller) != 0)
		err(1, "staging_benchmark: sched_getaffinity");

#pragma omp parallel num_threads(W)
	{
		int tid = omp_get_thread_num();
		int cpu = params.place.thread_cpu[1 + params.dicts_per_node + tid];
		if (params.bind_threads && pin_to_cpu(cpu) < 0)
			warn("staging_benchmark: cannot pin thread %d to CPU %d", tid, cpu);
		std::vector<u64> wc(n_dest * STAGING_LINE);       /* this thread's lines, its own first touch */
		std::vector<u8> fill(n_dest, 0);
		std::vector<u64> mine(cap * WORDS);               /* and its private buffer, for the floor */
		u64 seed = 0x9e3779b97f4a7c15ull * (tid + 1);

		/* every strategy starts from cursors at zero: a run reserved from anywhere else can straddle */
#pragma omp for schedule(static)
		for (u64 d = 0; d < n_dest; d++) {
			cur[d].next.store(0, std::memory_order_relaxed);
			for (u64 k = 0; k < cap * WORDS; k++)
				buf[d * cap * WORDS + k] = 0;
		}

		if (tid == 0)
			t0 = wtime();
#pragma omp barrier
		check[tid] ^= stage_local<WORDS>(mine.data(), cap, STAGING_POINTS, seed);
#pragma omp barrier
		if (tid == 0)
			rate[0] = W * STAGING_POINTS / (wtime() - t0);

#pragma omp for schedule(static)
		for (u64 d = 0; d < n_dest; d++)
			cur[d].next.store(0, std::memory_order_relaxed);
		if (tid == 0)
			t0 = wtime();
#pragma omp barrier
		check[tid] ^= stage_direct<WORDS>(buf.get(), cur.get(), dmask, cap, STAGING_POINTS, seed);
#pragma omp barrier
		if (tid == 0)
			rate[1] = W * STAGING_POINTS / (wtime() - t0);

#pragma omp for schedule(static)
		for (u64 d = 0; d < n_dest; d++)
			cur[d].next.store(0, std::memory_order_relaxed);
		if (tid == 0)
			t0 = wtime();
#pragma omp barrier
		check[tid] ^= stage_line<WORDS>(buf.get(), cur.get(), dmask, cap, STAGING_POINTS, seed,
		                                wc.data(), fill.data());
#pragma omp barrier
		if (tid == 0)
			rate[2] = W * STAGING_POINTS / (wtime() - t0);
	}

	if (sched_setaffinity(0, sizeof(caller), &caller) != 0)
		err(1, "staging_benchmark: sched_setaffinity");

	char hl[8], hd[8], hw[8];
	human_format(rate[0], hl);
	human_format(rate[1], hd);
	human_format(rate[2], hw);
	printf("    %8" PRIu64 " %8.1f MB %8.1f kB | %8s/s %8s/s %8s/s   %6.0f ns %6.0f ns\n", n_dest,
	       n_dest * cap * WORDS * 8 / 1e6, n_dest * STAGING_LINE * 8 / 1e3, hl, hd, hw,
	       1e9 * W / rate[1], 1e9 * W / rate[2]);
	for (int t = 0; t < W; t++)
		checksum ^= check[t];
}

/*
 * How fast can the producers of a node stage their own distinguished points, one shared buffer per
 * destination shard?  Prices the producer-side routing of PROBLEM.md §5 before it is built, at the
 * fan-out it will face -- `n_nodes * dicts_per_node` destinations, which is why the sweep is over
 * that.  Two different numbers come out of a row: the *node* rate is what replaces one comm thread's
 * 4.2 M DP/s, and the last two columns, what one thread spends per point, are the tax on a producer --
 * read them against the ~670 ns it spends walking a trail to produce the point (PROBLEM.md §9).
 *
 * What it does NOT measure: sealing a full buffer and handing it to the funnel.  The buffers here are
 * circular and are never sent, so a reservation that reaches the end simply wraps; the cache
 * footprint is honest, the slow path is missing.  The three strategies are timed in a fixed order on
 * warm buffers, so read `direct` against `line` -- they are adjacent -- rather than either against the
 * floor.  Rank-local, no MPI, no dictionary, no walking.
 */
template<typename Problem>
void staging_benchmark(const Problem& pb, const Options &opts, u64 nbytes_memory)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	Parameters params(opts, nbytes_memory, pb.n, pb.m);
	if (rank != 0)
		return;                     /* every rank is doing the same thing; one of them reports */

	u64 cap = 8;
	while (cap < params.buffer_capacity)
		cap *= 2;
	printf("Benchmarking DP staging: %d producer thread(s), %" PRIu64 " points per destination buffer,"
	       " %" PRIu64 " points per thread and measurement\n",
	       params.producers_per_node, cap, STAGING_POINTS);
	printf("  the sweep stops where the buffers no longer fit the %" PRIu64 " MB budget;"
	       " this job's own fan-out is %d (%d nodes x %d shards)\n",
	       nbytes_memory / 1000000, params.n_dicts, params.n_nodes, params.dicts_per_node);

	u64 checksum = 0;
	for (int words = 3; words >= 2; words--) {
		printf("  %d words per point (%d B%s), %" PRIu64 " points per write-combining line\n",
		       words, 8 * words, (words == 3) ? ", a length of its own" : ", the engine's DP",
		       STAGING_LINE / words);
		printf("        dest     buffer      line |      local      direct        line |"
		       " per point, per thread\n");
		for (u64 d = 1; ; d *= 4) {
			u64 shared_bytes = d * cap * words * 8;
			u64 private_bytes = (u64) params.producers_per_node * d * STAGING_LINE * 8;
			if (shared_bytes + private_bytes > nbytes_memory)
				break;
			if (words == 3)
				staging_row<3>(params, d, cap, checksum);
			else
				staging_row<2>(params, d, cap, checksum);
		}
	}
	printf("  checksum %016" PRIx64 "\n", checksum);
}
}
#endif
