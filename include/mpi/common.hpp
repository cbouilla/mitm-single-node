#ifndef MITM_MPI_COMMON
#define MITM_MPI_COMMON

#include <mpi.h>
#include <err.h>

#include "../common.hpp"

namespace mitm {

// enum tags {TAG_INTERCOMM, TAG_POINTS, TAG_SENDER_CALLHOME, TAG_RECEIVER_CALLHOME, TAG_ASSIGNMENT, TAG_SOLUTION};
// enum role {CONTROLLER, SENDER, RECEIVER, UNDECIDED};
// enum assignment {KEEP_GOING, NEW_VERSION};


class MpiParameters : public Parameters {
public:
	double ping_delay = 0.1;
	MPI_Comm comm;
	int mpi_rank, mpi_size;

	size_t buffer_capacity = 4096;               // somewhat arbitrary
	int n_recv_buffers = 32;                     // somewhat arbitrary
	int n_slack_buffers = 256;                   // somewhat arbitrary
    int min_backoff = 1000;                      // µs
    int max_backoff = 100000;                    // µs
    int max_current_send = 16;                   // somewhat arbitrary

	MpiParameters(MPI_Comm comm) : comm(comm)
	{
		MPI_Comm_rank(comm, &mpi_rank);
		MPI_Comm_size(comm, &mpi_size);
		verbose = (mpi_rank == 0);
	}
};

void BCast_result(const MpiParameters &params, vector<pair<u64,u64>> &result)
{
	// deal with the results
    vector<int> recvcounts(params.mpi_size);
    recvcounts[params.mpi_rank] = 2 * result.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, params.comm);

    vector<int> displs(params.mpi_size);
    int acc = 0;
    for (int i = 0; i < params.mpi_size; i++) {
        displs[i] = acc;
        acc += recvcounts[i];
    }

    vector<u64> tmp(acc);
    for (size_t i = 0; i < result.size(); i++) {
        int offset = displs[params.mpi_rank] + 2*i;
        auto [x0, x1] = result[i];
        tmp[offset] = x0;
        tmp[offset + 1] = x1;
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT, tmp.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, params.comm);
    result.clear();
    for (int i = 0; i < acc; i += 2)
        result.push_back(pair(tmp[i], tmp[i+1]));
}


#if 0
static void display_stats(u64 N, double start, int vlen, const MpiParameters &params)
{
	double rate = vlen * N / (wtime() - start);
    double rate_min = rate;
    double rate_max = rate;
    double rate_avg = rate;
    MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
    rate_avg /= params.size;
    double rate_std = (rate - rate_avg) * (rate - rate_avg);
    MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, params.local_comm);
    rate_std /= params.size;
    rate_std = std::sqrt(rate_std);
    if (params.rank == 0) {
        char hmin[8], hmax[8], havg[8], hstd[8];
        human_format(rate_min, hmin);
        human_format(rate_max, hmax);
        human_format(rate_avg, havg);
        human_format(rate_std, hstd);
        printf("Benchmark. f/s: min %s max %s avg %s std %s\n", hmin, hmax, havg, hstd);
    }
}

/* try to iterate for 1s. Return #it/s */
template<typename Problem>
void benchmark(const Problem& pb, const MpiParameters &params)
{
    if (params.rank == 0)
        printf("Benchmarking scalar implementation (using %d processes)\n", params.size);

    MPI_Barrier(params.world_comm);

    u64 N = 1ull << 26; 
    double start = wtime();
    u64 count = 0;
    for (u64 x = 0; x < N; x++) {
        u64 z = (x & 1) ? pb.f(x) : pb.g(x);
        u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
        int target = ((int) hash) % params.n_recv;
        if (target == 0)
            count += 1;
    }
    display_stats(N, start, 1, params);

    constexpr int vlen = Problem::vlen;
    if (vlen > 1) {
        if (params.rank == 0)
            printf("Benchmarking vector implementation (vlen=%d)\n", vlen);

        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen))); 
        u64 z[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choice[vlen];
        for (int i = 0; i < vlen; i++) {
        	choice[i] = i & 1;
            x[i] = i;
        }

		MPI_Barrier(params.world_comm);

        double start = wtime();
        u64 mask = make_mask(pb.n);
        u64 N = 1ull << 20; 
        for (u64 i = 0; i < N; i++) {
            pb.vfg(x, choice, z);
            for (int j = 0; j < vlen; j++)
                x[j] = z[j] & mask;
        }
        display_stats(N, start, vlen, params);
    }
}
#endif

}
#endif