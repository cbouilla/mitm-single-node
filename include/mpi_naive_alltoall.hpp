#ifndef MITM_MPI_NAIVE_ALLTOALL
#define MITM_MPI_NAIVE_ALLTOALL

#include <optional>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cmath>

#include "common.hpp"
#include "AbstractCollisionProblem.hpp"

#include <mpi.h>

/*
 * naive MITM w/ distributed dictionnary.  round-based Alltoallv version.
 */

namespace mitm {

template <bool EXPENSIVE_F, class AbstractProblem>
std::vector<std::pair<u64, u64>> naive_mpi_claw_search_alltoall(AbstractProblem &Pb, MpiParameters &params)
{
	static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
		"problem not derived from mitm::AbstractClawProblem");

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double start = wtime();
	u64 N = 1ull << Pb.n;
	CompactDict dict((1.25 * N) / size);
	std::vector<std::pair<u64, u64>> result;

    // expected #values received in each round by each process
	const u64 alpha = params.buffer_capacity;
	const u64 K = alpha*size;   // total values generated in each round
	// entries for each target follows binomial law (#values, 1/size);

	// generate K values :
	// each bucket <= alpha + sqrt(210 * alpha)   with proba 2^{-100}
	int limit = alpha + std::sqrt(120 * alpha);
	// if (rank == 0)
	//     printf("limit=%d\n", limit);

	u64 probe_false_pos = 0;
	u64 keys[3 * Pb.n];
	std::vector<u64> sendbuffer(size * limit);
	std::vector<u64> recvbuffer(size * limit);
	std::vector<int> sendcounts(size);
	std::vector<int> recvcounts(size);
	std::vector<int> displs(size);
	for (int i = 0; i < size; i++)
		displs[i] = limit * i;

	if (params.verbose) {
		char hbsize[8], hdsize[8];
		u64 bsize_process = 2 * sizeof(u64) * size * limit;
		u64 rank_per_node = size / params.n_nodes;
		human_format(bsize_process * rank_per_node, hbsize);
		u64 dsize_process = dict.n_slots * (sizeof(u64) + sizeof(u32));
		human_format(dsize_process * rank_per_node, hdsize);
		printf("RAM per node == %sB buffer + %sB dict\n", hbsize, hdsize);
	}

	for (int phase = 0; phase < 2; phase++) {
		// phase 0 == fill the dict with f()
		// phase 1 == probe the dict with g()
		if (rank == 0)
			printf("Starting phase %d\n", phase);
		double phase_start = wtime();
		double last_display = phase_start;
		double wait = 0;

		const u64 nrounds = (N + K - 1) / K;
		for (u64 round = 0; round < nrounds; round ++) {
			// reset sendcounts
			for (int i = 0; i < size; i++)
				sendcounts[i] = 0;

			// round = [N * round / nrounds : N * (round + 1) / nrounds]
			u64 round_lo = N * round / nrounds;
			u64 round_hi = N * (round + 1) / nrounds;
			u64 round_size = round_hi - round_lo;
			u64 process_lo = round_lo + rank * round_size / size;
			u64 process_hi = round_lo + (rank + 1) * round_size / size;
			for (u64 x = process_lo; x < process_hi; x++) {
				u64 z = (phase == 0) ? Pb.f(x) : Pb.g(x);
				u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
				int target = ((int) hash) % size;
				assert(sendcounts[target] < limit);
				u64 offset = limit * target + sendcounts[target];
				if (EXPENSIVE_F) {
					sendbuffer[offset] = x;
					sendbuffer[offset + 1] = z;
					sendcounts[target] += 2;
				} else {
					sendbuffer[offset] = x;
					sendcounts[target] += 1;
				}
			}
			double start_comm = wtime();

			// exchange buffer sizes;
			MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

			// exchange data / TODO: IAlltoallv
			MPI_Alltoallv(sendbuffer.data(), sendcounts.data(), displs.data(), MPI_UINT64_T, 
					  recvbuffer.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, MPI_COMM_WORLD);

			wait += wtime() - start_comm;

			for (int i = 0; i < size; i++)
				for (int j = 0; j < recvcounts[i]; j++) {
					u64 x = recvbuffer[i * limit + j];
					u64 z;
					if (EXPENSIVE_F) {
						j += 1;
						z = recvbuffer[i * limit + j];
					} else {
						z = (phase == 0) ? Pb.f(x) : Pb.g(x);
					}
					if (phase == 0) {
						// inject data into dict
						dict.insert(z, x);
					} else {
						// probe dict
						int nkeys = dict.probe(z, keys);
						for (int k = 0; k < nkeys; k++) {
							u64 y = keys[k];
							if (z != Pb.f(y)) {
								probe_false_pos += 1;
								continue;    // false positive from truncation in the hash table
							}
							if (Pb.is_good_pair(y, x)) {
								// printf("\nfound golden collision !!!\n");
								result.push_back(std::pair(y, x));
							}
						}
					}
				}

			// verbosity
			double now = wtime();
			if (rank == 0 && now - last_display >= 0.5) {
				char frate[8], nrate[8];
				double delta = now - phase_start;
				last_display = now;
				human_format(K * (round + 1) / delta, frate);
				human_format(8 * K * (round + 1) / delta, nrate);
				printf("Round %" PRId64 " / %" PRId64 ".  Wait/round = %.3fs (%.1f%%).  %s f()/s.  Net=%sB/s\n",
				   round, nrounds, wait/(1+round), 100.*wait / delta, frate, nrate);
				fflush(stdout);
			}
		} //  round
		if (rank == 0)
			printf("Phase: %.1fs\n", wtime() - phase_start);
	} // phase
	if (rank == 0)
		printf("Total: %.1fs\n", wtime() - start);
	BCast_result(params, result);
	return result;
}

}

#endif
