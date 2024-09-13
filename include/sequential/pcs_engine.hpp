#ifndef MITM_ENGINE_SEQ
#define MITM_ENGINE_SEQ

#include <cmath>
#include <cstdio>

#include "common.hpp"
#include "engine_common.hpp"

namespace mitm {

/* The sequence of mixing function and of evaluation points is deterministic given `prng` */

class ScalarSequentialEngine : Engine {
public:

template<class ProblemWrapper>
static optional<tuple<u64,u64,u64>> run(ProblemWrapper& wrapper, Parameters &params, PRNG &prng)
{
    int jbits = std::log2(10 * params.w) + 8;
    u64 jmask = make_mask(jbits);
    u64 w = PcsDict::get_nslots(params.nbytes_memory, 1);
    PcsDict dict(jbits, w);
    
    Counters ctr;
    ctr.ready(wrapper.n, w);

    double log2_w = std::log2(w);
    printf("Starting collision search with seed=%016" PRIx64 " (scalar engine)\n", prng.seed);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, params.points_per_version, std::log2(params.points_per_version));

    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    for (u64 nver = 0; nver < params.max_versions; nver++) {
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
            
            solution = process_distinguished_point(wrapper, ctr, params, dict, i, root_seed, j, end, len);
            if (solution)
                break;
        }
        dict.flush();
        ctr.flush_dict();
        if (solution)
            break;
    }
    ctr.done();
    return solution;
}
};

class VectorSequentialEngine : Engine {
public:

static void start_chain(const Parameters &params, u64 out_mask, u64 root_seed, u64 &j, u64 x[], u64 len[], u64 seed[], int k)
{
    u64 start;
    for (;;) {
        j += 1;
        start = (root_seed + j * params.multiplier) & out_mask;
        if (not is_distinguished_point(start, params.threshold))  // refuse to start from a DP
            break;
    }
    x[k] = start;
    len[k] = 0;
    seed[k] = j;    
}


template<class ProblemWrapper>
static optional<tuple<u64,u64,u64>> run(ProblemWrapper& wrapper, Parameters &params, PRNG &prng)
{
    int jbits = std::log2(10 * params.w) + 8;
    u64 w = PcsDict::get_nslots(params.nbytes_memory, 1);
    PcsDict dict(jbits, w);

    Counters ctr;
    ctr.ready(wrapper.n, w);

    double log2_w = std::log2(w);
    printf("Starting collision search with seed=%016" PRIx64 " (vectorized engine)\n", prng.seed);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, params.points_per_version, std::log2(params.points_per_version));

    optional<tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    u64 i, root_seed, j;
    ctr.n_dp_i = params.points_per_version;   // trigger new version right from the start
    constexpr int vlen = ProblemWrapper::vlen;
    u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
    u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
    u64 len[vlen], seed[vlen];

    for (;;) {
        /* These simulations show that if 10w distinguished points are generated
         * for each version of the function, and theta = 2.25sqrt(w/n) then ...
         */
        if (ctr.n_dp_i >= params.points_per_version) {
            /* new version of the function */
            i = prng.rand() & wrapper.out_mask;
            root_seed = prng.rand();
            j = 0;
            dict.flush();
            ctr.flush_dict();
            /* restart all the chains */
            for (int k = 0; k < vlen; k++)
                start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, k);
        }

        /* advance all the chains */
        wrapper.vmixf(i, x, y);

        /* test for distinguished points */ 
        for (int k = 0; k < vlen; k++) {
            len[k] += 1;
            x[k] = y[k];
            bool dp = is_distinguished_point(x[k], params.threshold);
            bool failure = (not dp && len[k] == params.dp_max_it);
            if (failure)
                ctr.dp_failure();
            if (dp) {
                ctr.found_distinguished_point(len[k]);
                auto solution = process_distinguished_point(wrapper, ctr, params, dict, i, root_seed, seed[k], x[k], len[k]);
                if (solution)
                    return *solution;
            }
            if (dp || failure)
                start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, k);
        }
    } // main loop
    return nullopt;
}
};


}
#endif
