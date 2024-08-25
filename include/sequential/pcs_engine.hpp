#ifndef MITM_ENGINE_SEQ
#define MITM_ENGINE_SEQ

#include <cmath>
#include <cstdio>

#include "common.hpp"
#include "engine_common.hpp"

namespace mitm {

class SequentialEngine : Engine {
public:

/* The sequence of mixing function and of evaluation points is deterministic given `prng` */
template<class ProblemWrapper>
static tuple<u64,u64,u64> run(ProblemWrapper& wrapper, Parameters &params, PRNG &prng)
{
    int jbits = std::log2(10 * params.w) + std::log2(1 / params.theta) + 8;
    u64 jmask = make_mask(jbits);
    u64 w = PcsDict::get_nslots(params.nbytes_memory, 1);
    PcsDict dict(jbits, w);
    
    Counters ctr;
    ctr.ready(wrapper.n, w);

    double log2_w = std::log2(w);
    printf("Starting collision search with seed=%016" PRIx64 "\n", prng.seed);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, params.points_per_version, std::log2(params.points_per_version));

    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    for (;;) {
        /* These simulations show that if 10w distinguished points are generated
         * for each version of the function, and theta = 2.25sqrt(w/n) then ...
         */
        u64 i = prng.rand() & wrapper.out_mask;           /* index of families of mixing functions */
        u64 root_seed = prng.rand();
        u64 j = 0;
        while (ctr.n_dp_i < params.points_per_version) {
            j += 1;
            
            assert((j & jmask) == j);

            /* start a new chain from a fresh "random" non-distinguished starting point */
            u64 start = (root_seed + j * params.multiplier) & wrapper.out_mask;
            if (is_distinguished_point(start, params.threshold))  // refuse to start from a DP
                continue;

            auto dp = generate_dist_point(wrapper, i, params, start);
            if (not dp) {
                ctr.dp_failure();
                continue;
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);
            
            auto solution = process_distinguished_point(wrapper, ctr, params, dict, i, root_seed, j, end, len);
            if (solution)
                return *solution;
        }
        dict.flush();
        ctr.flush_dict();
    } // main loop
}
};

}
#endif
