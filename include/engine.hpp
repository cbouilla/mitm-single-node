#ifndef MITM_ENGINE
#define MITM_ENGINE

#include <cmath>
#include <cassert>
#include <cstdio>
#include <optional>
#include <tuple>

#include <omp.h>

#include "AbstractCollisionProblem.hpp"
#include "common.hpp"
#include "counters.hpp"
#include "dict.hpp"

namespace mitm {

class Parameters {
public:
    double difficulty = -1;   /* -log2(proportion of distinguished points). -1 == auto-choose */
    u64 nbytes_memory = 0;    /* how much RAM for the dictionary. 0 == use everything */
    double beta = 10;         /* use function variant for beta*w distinguished points */
    u64 threshold;            /* any integer less than this is a DP */
    u64 dp_max_it;            /* how many iterations to find a DP */

    static double optimal_theta(double log2_w, int n)
    {
        return std::log2(2.25) + 0.5 * (log2_w - n);
    }

    void finalize(int n)
    {
         /* If the user did not specify how much memory to be use, use all the RAM */
        if (nbytes_memory == 0) {
            nbytes_memory = get_available_memory();
            char hm[8];
            human_format(nbytes_memory, hm);
            printf("AUTO-TUNING: using %sB of RAM\n", hm);
        }
        u64 nslots = Dict<std::pair<u64, u64>>::get_nslots(nbytes_memory);
        double log2_w = std::log2(nslots);
        double theta = optimal_theta(log2_w, n);
        /* auto-choose the difficulty if not set */
        if (difficulty < 0) {
            difficulty = -theta;
            if (difficulty < 0)
                difficulty = 0;
            printf("AUTO-TUNING: setting difficulty %.2f\n", difficulty);
        } else {
            printf("NOTICE: using difficulty=%.2f vs ``optimal''=%.2f\n", difficulty, theta);
        }
        threshold = std::pow(2., n - difficulty);
        dp_max_it = 20 * std::pow(2., difficulty);

        /* display warnings if problematic choices were made */
        if (difficulty < 0) {
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
            printf("---> negative difficulty (reset to 0)\n");
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");            
        }
        if (n - difficulty < log2_w) {
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
            printf("---> Too much memory (2^%.2f slots) but only 2^%.2f possible distinguished points\n",
                log2_w, n - difficulty);
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
        }
    }
};


/******************************************************************************/

inline bool is_distinguished_point(u64 x, u64 threshold)
{
    return x <= threshold;
}


/*
 * Given an input, iterate functions either F or G until a distinguished point
 * is found, save the distinguished point in out_pt and output_bytes
 * Then return `true`. If the iterations limit is passed, returns `false`.
 */
template<typename ConcreteProblem>
std::optional<std::pair<u64,u64>> generate_dist_point(const ConcreteProblem& Pb, u64 i, const Parameters &params, u64 x)
{
    /* The probability, p, of NOT finding a distinguished point after the loop is
     * Let: theta := 2^-d
     * difficulty, N = k*2^difficulty then,
     * p = (1 - theta)^N =>  let ln(p) <= -k
     */
    for (u64 j = 0; j < params.dp_max_it; j++) {
        u64 y = Pb.mixf(i, x);
        if (is_distinguished_point(y, params.threshold))
            return std::make_optional(std::pair(y, j+1));
        x = y;
    }
    return std::nullopt; /* no distinguished point was found after too many iterations */
}



/* Given two inputs that lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * add a drawing to illustrate this.
 */
template<typename ConcreteProblem>
std::optional<std::tuple<u64,u64,u64>> walk(const ConcreteProblem& Pb, u64 i, u64 x0, u64 len0, u64 x1, u64 len1)
{
  /****************************************************************************+
   *            walk the longest sequence until they are equal                 |
   * Two chains that leads to the same distinguished point but not necessarily |
   * have the same length. e.g.                                                |
   *                                                                           |
   * chain1: ----------------x-------o                                         |
   *                        /                                                  |
   *          chain2: ------                                                   |
   *                                                                           |
   * o: is a distinguished point                                               |
   * x: the collision we're looking for                                        |
   ****************************************************************************/

  /* move the longest sequence until the remaining number of steps is equal */
  /* to the shortest sequence. */
  for (; len0 > len1; len0--)
    x0 = Pb.mixf(i, x0);
  for (; len0 < len1; len1--)
    x1 = Pb.mixf(i, x1);
  
  /* now both sequences needs exactly `len` steps to reach distinguished point */
  for (u64 j = 0; j < len0; ++j) {
    /* walk them together and check each time if their output are equal     */
    /* return as soon equality is found. The equality could be a robinhood. */
    u64 y0 = Pb.mixf(i, x0);
    u64 y1 = Pb.mixf(i, x1);

    /* First, do the outputs collide? If yes, return true and exit. */
    if (y0 == y1) {
      /* careful: x0 & x1 contain inputs before mixing */
      return std::make_optional(std::tuple(x0, x1, y0));
    }
    x0 = y0;
    x1 = y1;
  }
  return std::nullopt; /* we did not find a common point */
}

template<typename ConcreteProblem>
double benchmark(const ConcreteProblem& Pb)
{
    u64 x = 0;
    u64 i = 42;
    u64 threshold = (1ull << (Pb.n - 1));
    double start = wtime();
    u64 N = 1000000;
    u64 total = 0;
    #pragma omp parallel reduction(+:total)
    {
        total = N;
        for (u64 j = 0; j < N; j++) {
            u64 y = Pb.mixf(i, x);
            if (is_distinguished_point(y, threshold))
                x = j;
            else
                x = y;
        }
    }
    return total / (wtime() - start);
}


/*
 *  The sequence of mixing function is deterministic given `prng`  
 */
template<typename ConcreteProblem>
std::tuple<u64,u64,u64> search_generic(const ConcreteProblem& Pb, Parameters &params, PRNG &prng)
{
    Counters &ctr = Pb.ctr; 
    params.finalize(Pb.n);

    printf("Benchmarking... ");
    fflush(stdout);
    double it_per_s = benchmark(Pb);
    char hitps[8];
    human_format(it_per_s, hitps);
    printf("%s iteration/s (using all cores)\n", hitps);

    printf("Starting collision search with seed=%016" PRIx64 ", difficulty=%.2f\n", prng.seed, params.difficulty);

    Dict<std::pair<u64, u64>> dict(params.nbytes_memory);
    double log2_w = std::log2(dict.n_slots);
    printf("Initialized a dict with %" PRId64 " slots = 2^%0.2f slots\n", dict.n_slots, log2_w);
    
    ctr.ready(Pb.n, dict.n_slots);

    printf("Expected iterations / collision = (2^%0.2f + 2^%.2f) \n", 
        Pb.n - params.difficulty - log2_w, 1 + params.difficulty);

    printf("Expected #iterations = (2^%0.2f + 2^%.2f) \n", 
        (Pb.n - 1) + (Pb.n - params.difficulty - log2_w), Pb.n + params.difficulty);


    u64 i = 0;                 /* index of families of mixing functions */
    std::optional<std::tuple<u64,u64,u64>> solution;    /* (i, x0, x1)  */
    const u64 points_per_version = params.beta * dict.n_slots;
    u64 n_dist_points = 0;     /* #DP found with this i */

    printf("Generating %.1f*w = %" PRId64 " = 2^%0.2f distinguished point / version\n", 
        params.beta, points_per_version, std::log2(points_per_version));

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
            if (tmp > points_per_version)
                break;

            /* start a new chain from a fresh "random" starting point */
            u64 start0 = (root_seed + n_dist_points) & Pb.mask;
            auto dp = generate_dist_point(Pb, i, params, start0);
            if (not dp) {
                ctr.dp_failure();
                continue;       /* bad chain start */
            }

            #pragma omp atomic
            n_dist_points += 1;

            auto [end, len0] = *dp;
            ctr.found_distinguished_point(len0);
            auto probe = dict.pop_insert(end, std::pair(start0, len0));
            if (not probe) {
                ctr.probe_failure();
                continue;
            }

            auto [start1, len1] = *probe;
            auto collision = walk(Pb, i, start0, len0, start1, len1);  
            if (not collision) {
                ctr.walk_failure();
                continue;         /* robin-hood, or dict false positive */
            }

            auto [x0, x1, y] = *collision;
            if (x0 == x1) {
                ctr.collision_failure();
                continue;    /* duh */
            }

            ctr.found_collision();
            if (Pb.mix_good_pair(i, x0, x1)) {
                printf("\nFound golden collision!\n");
                solution = std::make_optional(std::tuple(i, x0, x1));
            }
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

}
#endif
