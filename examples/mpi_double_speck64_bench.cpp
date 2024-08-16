#include <cassert>
#include <getopt.h>
#include <err.h>

#include <mpi.h>

#include "mpi_common.hpp"
#include "mpi_naive_alltoall.hpp"
#include "mpi_naive_isend.hpp"
#include "double_speck64_problem.hpp"

u64 seed = 1337;

namespace mitm {

/* try to iterate for 1s. Return #it/s */
template<typename AbstractProblem>
static void naive_benchmark(const AbstractProblem& Pb, MpiParameters &params)
{
    MPI_Barrier(params.world_comm);

    u64 N = 1ull << 30; 
    double start = wtime();
    u64 count = 0;
    for (u64 x = 0; x < N; x++) {
        u64 z = (x & 1) ? Pb.f(x) : Pb.g(x);
        u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
        int target = ((int) hash) % params.n_recv;
        if (target == 0)
            count += 1;
    }
    double rate = N / (wtime() - start);
    double rate_min = rate;
    double rate_max = rate;
    double rate_avg = rate;
    MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
    rate_avg /= params.size;
    double rate_std = (rate - rate_avg) * (rate - rate_avg);
    MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, params.local_comm);
    rate_std = std::sqrt(rate_std);
    if (params.rank == 0) {
        char hmin[8], hmax[8], havg[8], hstd[8];
        human_format(rate_min, hmin);
        human_format(rate_max, hmax);
        human_format(rate_avg, havg);
        human_format(rate_std, hstd);
        printf("Benchmark. f/s: min %s max %s avg %s std %s (aggregated over %d processes)\n", hmin, hmax, havg, hstd, params.size);
    }
}
}

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    
    mitm::MpiParameters params;
    params.setup(MPI_COMM_WORLD, 0);  // no controller process
    
    mitm::PRNG prng(seed);
    mitm::DoubleSpeck64_Problem Pb(32, prng);
    mitm::naive_benchmark(Pb, params);

    MPI_Finalize();    
    return EXIT_SUCCESS;
}
