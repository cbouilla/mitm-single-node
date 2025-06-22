#ifndef MITM_MPI_BENCHMARK
#define MITM_MPI_BENCHMARK

#include <mpi.h>
#include <err.h>

#include "../common.hpp"

namespace mitm {

static void display_stats(u64 N, double start, int vlen, const MpiParameters &params)
{
	double rate = vlen * N / (wtime() - start);
    double rate_min = rate;
    double rate_max = rate;
    double rate_avg = rate;
    MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, params.comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, params.comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, params.comm);
    rate_avg /= params.mpi_size;
    double rate_std = (rate - rate_avg) * (rate - rate_avg);
    MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, params.comm);
    rate_std /= params.mpi_size;
    rate_std = std::sqrt(rate_std);
    if (params.verbose) {
        char hmin[8], hmax[8], havg[8], hstd[8];
        human_format(rate_min, hmin);
        human_format(rate_max, hmax);
        human_format(rate_avg, havg);
        human_format(rate_std, hstd);
        printf("Benchmark. f/s (per host): min %s max %s avg %s std %s\n", hmin, hmax, havg, hstd);
    }
}

/* try to iterate for 1s. Return #it/s */
template<typename Problem>
void benchmark(const Problem& pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    assert(pb.n >= 26);

    if (params.verbose)
        println("Claw-finding: {{0,1}}^{} --> {{0,1}}^{}", pb.n, pb.m);

    if (params.n_threads <= 0) {
        params.n_threads = omp_get_max_threads();
        if (params.verbose)
            println("Autodetect: using {} threads", params.n_threads);
    }
    if (params.verbose)
        println("Benchmarking scalar implementation (using {} MPI processes with {} threads)", params.mpi_size, params.n_threads);

    MPI_Barrier(params.comm);

    u64 N = 1ull << 26; 
    double start = wtime();
    #pragma omp parallel
    {
        u64 count = 0;
        for (u64 x = 0; x < N; x++) {
            u64 z = (x & 1) ? pb.f(x) : pb.g(x);
            u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
            int target = ((int) hash) % params.mpi_size;
            if (target == 0)
                count += 1;
        }
    }
    display_stats(N, start, 1, params);

    constexpr int vlen = Problem::vlen;
    if (vlen > 1) {
        if (params.verbose)
            println("Benchmarking vector implementation (vlen={})", vlen);

        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen))); 
        u64 z[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choice[vlen];
        for (int i = 0; i < vlen; i++) {
        	choice[i] = i & 1;
            x[i] = i;
        }

		MPI_Barrier(params.comm);
        
        start = wtime();
        u64 N = 1ull << 20; 
        #pragma omp parallel
        {
            u64 mask = make_mask(pb.n);
            for (u64 i = 0; i < N; i++) {
                pb.vfg(x, choice, z);
                for (int j = 0; j < vlen; j++)
                    x[j] = z[j] & mask;
            }
        }
        display_stats(N, start, vlen, params);
    }
}

}
#endif