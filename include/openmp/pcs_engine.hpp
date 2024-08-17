#ifndef MITM_ENGINE_OMP
#define MITM_ENGINE_OMP

#include <cmath>
#include <cstdio>

#include <omp.h>

#include "common.hpp"
#include "engine_common.hpp"


/* multi-threaded collision-finding engine using OpenMP */

namespace mitm {

/* Non-essential counters but helpful to have, e.g. n_collisions/sec */
class OmpCounters : public BaseCounters {
public:
    OmpCounters() {}
    OmpCounters(bool display_active) : BaseCounters(display_active) {}

    void dp_failure()
    { 
        #pragma omp atomic
        bad_dp += 1; 
    }
    
    void probe_failure()
    {
        #pragma omp atomic      
        bad_probe += 1;
    }
    
    void walk_failure() {
        #pragma omp atomic      
        bad_walk += 1;
    }
    
    void collision_failure() {
        #pragma omp atomic      
        bad_collision += 1;
    }

    void function_evaluation()
    {
        #pragma omp atomic                // PROBLEM! This likely entails a significant performance penalty...
        n_points += 1;
    }

    // call this when a new DP is found
    void found_distinguished_point(u64 chain_len)
    {
        #pragma omp atomic
        n_dp += 1;
        #pragma omp atomic
        n_dp_i += 1;
        #pragma omp atomic
        n_points_trails += chain_len;
        if ((n_dp % interval == interval - 1))
             #pragma omp critical(counters)
             display();
    }

    void found_collision() {
        #pragma omp atomic
        n_collisions += 1;
        #pragma omp atomic
        n_collisions_i += 1;
    }
};  


class OpenMPEngine {
public:
using Counters = OmpCounters;

/* try to iterate for 1s. Return #it/s */
template<typename ConcreteProblem>
static double benchmark(const ConcreteProblem& Pb)
{
    double rate = 0;
    #pragma omp parallel reduction(+:rate)
    {
        rate = sequential_benchmark(Pb);
        #pragma omp critical
        printf("thread %d, rate=%.0f it/s\n", omp_get_thread_num(), rate);
    }
    return rate;
}


/* The sequence of mixing function is deterministic given `prng` */
template<typename ConcreteProblem>
static tuple<u64,u64,u64> run(const ConcreteProblem& Pb, Parameters &params, PRNG &prng)
{
    Counters &ctr = Pb.ctr;

    printf("Benchmarking... ");
    fflush(stdout);
    double it_per_s_seq = sequential_benchmark(Pb);
    double it_per_s = benchmark(Pb);
    // NOT DRY wrt MPI
    char hitps[8], hitps_seq[8];
    human_format(it_per_s_seq, hitps_seq);
    human_format(it_per_s, hitps);

    int nthreads = omp_get_max_threads();
    printf("%s it/s (one thread).  %s it/s (%d threads)\n", hitps_seq, hitps, nthreads);

    PcsDict dict(params.nbytes_memory);
    params.finalize(Pb.n, dict.n_slots);
    ctr.ready(Pb.n, dict.n_slots);
    double log2_w = std::log2(dict.n_slots);

    printf("Starting collision search with seed=%016" PRIx64 ", difficulty=%.2f\n", prng.seed, params.difficulty);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    printf("Expected iterations / collision = (2^%0.2f + 2^%.2f) \n", 
        Pb.n - params.difficulty - log2_w, 1 + params.difficulty);
    printf("Expected #iterations = (2^%0.2f + 2^%.2f) \n", 
        (Pb.n - 1) + (Pb.n - params.difficulty - log2_w), Pb.n + params.difficulty);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, params.points_per_version, std::log2(params.points_per_version));

    u64 i = 0;                 /* index of families of mixing functions */
    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    u64 n_dist_points = 0;     /* #DP found with this i */
    u64 root_seed = prng.rand();

    #pragma omp parallel    
    for (;;) {
        /* These simulations show that if 10w distinguished points are generated
         * for each version of the function, and theta = 2.25sqrt(w/n) then ...
         */

        for (;;) {
            u64 tmp;
            #pragma omp atomic read
            tmp = n_dist_points;
            if (tmp > params.points_per_version || solution)
                break;

            #pragma omp atomic
            n_dist_points += 1;

            /* start a new chain from a fresh "random" starting point */
            u64 start = (root_seed + n_dist_points) & Pb.mask;
            auto dp = generate_dist_point(Pb, i, params, start);
            if (not dp) {
                ctr.dp_failure();
                continue;       /* bad chain start ------ FIXME infinite loop if this happen */
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);
            
            auto outcome = process_distinguished_point(Pb, ctr, dict, i, start, end, len);
            if (outcome)
                solution = outcome;
        }

        if (solution)
            break;

        #pragma omp barrier

        /* change the mixing function */
        #pragma omp single
        {
            n_dist_points = 0;
            i = prng.rand() & Pb.mask; 
            dict.flush();
            ctr.flush_dict();
            root_seed = prng.rand();
        }
    } // main loop
    return *solution;
}
};

}
#endif
