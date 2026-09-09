#ifndef MITM_BENCHMARK
#define MITM_BENCHMARK

#include <mpi.h>
#include <cmath>
#include <type_traits>

#include "tools.hpp"
#include "parameters.hpp"
#include "problem.hpp"

/*
 * How fast do a problem's functions iterate, per rank and over the ranks?  What every driver's
 * --benchmark runs; needs no RAM budget, only the Options.  The dictionary and staging benchmarks
 * that used to live here were built on PCS's parameters and come back when PCS does; router_bench
 * measures the communication path on its own.
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
	if (rank == 0)
		fmt::print("Benchmark. f/s: min {} max {} avg {} std {}\n", human_format(rate_min),
		           human_format(rate_max), human_format(rate_avg), human_format(rate_std));
}

/*
 * Iterate the scalar functions 2^26 times, then the vector ones 2^20 times if vlen > 1, on every rank,
 * and print the rates.  A claw problem alternates f and g, a collision problem iterates f alone.  Both
 * loops walk ONE dependent chain per lane, the shape a trail and a collision resolution both have, so
 * that the two rates can be divided: what that ratio is worth is how much a producer gains by batching
 * its resolutions.  Each loop's last value is printed, and that is the only reason it survives -O3:
 * nothing else here is observable.  Collective
 */
template<typename Problem>
void benchmark(const Problem& pb, const Options &opts)
{
	constexpr bool claw = std::is_base_of<AbstractClawProblem, Problem>::value;
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
	for (u64 i = 0; i < N; i++) {
		if constexpr (claw)
			x1 = ((i & 1) ? pb.f(x1) : pb.g(x1)) & mask;
		else
			x1 = pb.f(x1) & mask;
	}
	display_stats(N, start, 1, opts.mpi_comm, rank, n_nodes);
	if (rank == 0)
		printf("  checksum %016" PRIx64 "\n", x1);

	constexpr int vlen = Problem::vlen;
	if constexpr (vlen > 1) {
		if (rank == 0)
			printf("Benchmarking vector implementation (vlen=%d)\n", vlen);

		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 z[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		bool choice[vlen];      /* what vfg selects, for a claw problem: half the lanes f, the others g */
		for (int i = 0; i < vlen; i++) {
			choice[i] = i & 1;
			x[i] = i;
		}

		MPI_Barrier(opts.mpi_comm);

		double start = wtime();
		u64 N = 1ull << 20;
		for (u64 i = 0; i < N; i++) {
			if constexpr (claw)
				pb.vfg(x, choice, z);
			else
				pb.vf(x, z);
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
}
#endif
