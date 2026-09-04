#ifndef MITM_BENCHMARK
#define MITM_BENCHMARK

#include <mpi.h>
#include <cmath>

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
 * have, so that the two rates can be divided: what that ratio is worth is how much a walker gains by
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
 * How fast can an inserter probe its shard, with nothing else running?  R_p has only ever been
 * measured through the whole pipeline, tangled up with the queues and the router (PROBLEM.md §6).
 * This is the shard alone, at the size and the placement the engine would give it: the rank's
 * `inserters_per_node` threads, pinned exactly where Parameters puts them, each hammering its own
 * PcsDict with pseudo-random endpoints.  The rate is reported per quarter-`w` batch, because the
 * probe's cost depends on how full the shard is: a round starts empty and ends at `beta`
 * insertions per slot.  Rank-local, no MPI, no walkers, no comm thread.
 */
template<typename Problem>
void probe_benchmark(const Problem& pb, const Options &opts, u64 nbytes_memory)
{
	int rank;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	Parameters params(opts, nbytes_memory, pb.n, pb.m);
	if (rank != 0)
		return;                     /* every rank is doing the same thing; one of them reports */

	int I = params.inserters_per_node;
	char hw[8], hws[8];
	human_format(params.w, hw);
	human_format(params.w_shard, hws);
	printf("Benchmarking dictionary probes: %d shard(s) of %s slots (%s total, %.1f MB each)\n",
	       I, hws, hw, params.w_shard * sizeof(u64) / 1e6);

	u64 batch = params.w_shard / 4;                   /* probes per timed batch, per shard */
	int n_batches = std::max(1, (int) (4 * params.beta));   /* --beta bounds the run: it is long */
	u64 key_range = std::ldexp(1., pb.n) / params.n_inserters;
	std::vector<double> rate(n_batches * I);
	std::vector<u64> hits(I, 0);

#pragma omp parallel num_threads(I)
	{
		int tid = omp_get_thread_num();
		if (params.bind_threads && pin_to_cpu(params.thread_cpu[1 + tid]) < 0)
			warn("probe_benchmark: cannot pin thread %d to CPU %d", tid, params.thread_cpu[1 + tid]);
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
	u64 total_hits = 0;
	for (int t = 0; t < I; t++)
		total_hits += hits[t];
	printf("  %" PRIu64 " hits in %" PRIu64 " probes (%.2f%%)\n", total_hits,
	       (u64) n_batches * batch * I, 100. * total_hits / ((double) n_batches * batch * I));
}

}
#endif
