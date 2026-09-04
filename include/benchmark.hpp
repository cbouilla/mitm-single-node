#ifndef MITM_BENCHMARK
#define MITM_BENCHMARK

#include <mpi.h>
#include <cmath>

#include "tools.hpp"
#include "parameters.hpp"

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

/* iterate f / g 2^26 times, then vfg 2^20 times if vlen > 1, on every rank, and print the rates.  Collective */
template<typename Problem>
void benchmark(const Problem& pb, const Options &opts)
{
	int rank, n_nodes;
	MPI_Comm_rank(opts.mpi_comm, &rank);
	MPI_Comm_size(opts.mpi_comm, &n_nodes);
	int n_inserters = n_nodes * opts.inserters_per_node;

	if (rank == 0)
		printf("Benchmarking scalar implementation (using %d processes)\n", n_nodes);

	MPI_Barrier(opts.mpi_comm);

	u64 N = 1ull << 26;
	double start = wtime();
	u64 count = 0;
	for (u64 x = 0; x < N; x++) {
		u64 z = (x & 1) ? pb.f(x) : pb.g(x);
		u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
		int target = ((int) hash) % n_inserters;
		if (target == 0)
			count += 1;
	}
	display_stats(N, start, 1, opts.mpi_comm, rank, n_nodes);

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
		u64 mask = make_mask(pb.n);
		u64 N = 1ull << 20;
		for (u64 i = 0; i < N; i++) {
			pb.vfg(x, choice, z);
			for (int j = 0; j < vlen; j++)
				x[j] = z[j] & mask;
		}
		display_stats(N, start, vlen, opts.mpi_comm, rank, n_nodes);
	}
}

}
#endif
