#ifndef MITM_ENGINE_SEQ
#define MITM_ENGINE_SEQ

#include <cmath>
#include <cstdio>

#include "common.hpp"
#include "engine_common.hpp"

namespace mitm {

class SequentialEngine : Engine {
public:

/* The sequence of mixing function is deterministic given `prng` */
template<class ProblemWrapper>
static tuple<u64,u64,u64> run(ProblemWrapper& wrapper, Parameters &params, PRNG &prng)
{
    Counters ctr;
    int n = wrapper.pb.n;
    PcsDict dict(params.nbytes_memory);
    ctr.ready(n, dict.n_slots);
    double log2_w = std::log2(dict.n_slots);

    printf("Starting collision search with seed=%016" PRIx64 "\n", prng.seed);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, params.points_per_version, std::log2(params.points_per_version));

    u64 mask = make_mask(wrapper.pb.m);
    u64 i = 0;                 /* index of families of mixing functions */
    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    u64 root_seed = prng.rand();

    while (not solution) {
        /* These simulations show that if 10w distinguished points are generated
         * for each version of the function, and theta = 2.25sqrt(w/n) then ...
         */

        u64 j = 0;
        while (ctr.n_dp_i < params.points_per_version) {
            j += 0x2545f4914f6cdd1dull;
            
            /* start a new chain from a fresh "random" starting point */
            u64 start = (root_seed + j) & mask;
            if (is_distinguished_point(start, params.threshold))  // refuse to start from a DP
                continue;

            auto dp = generate_dist_point(wrapper, i, params, start);
            if (not dp) {
                ctr.dp_failure();
                continue;
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);
            
            auto outcome = process_distinguished_point(wrapper, ctr, dict, i, start, end, len);
            if (outcome) {
                solution = outcome;
                break;
            }
        }

        /* change the mixing function */
        i = prng.rand() & mask; 
        dict.flush();
        ctr.flush_dict();
        root_seed = prng.rand();

    } // main loop
    return *solution;
}
};

}
#endif
