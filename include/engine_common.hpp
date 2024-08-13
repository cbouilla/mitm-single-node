#ifndef MITM_ENGINE_COMMON
#define MITM_ENGINE_COMMON

#include <cmath>
#include <cassert>
#include <cstdio>
#include <optional>
#include <tuple>

#include "common.hpp"
#include "AbstractCollisionProblem.hpp"
#include "counters.hpp"
#include "dict.hpp"

namespace mitm {

class Parameters {
public:
    u64 nbytes_memory = 0;        /* how much RAM to use on each machine */
    double difficulty = -1;       /* -log2(proportion of distinguished points). -1 == auto-choose */
    double beta = 10;             /* use function variant for beta*w distinguished points */
    int n_nodes = 1;              /* #hosts (with shared RAM) */
    bool verbose = 1;             /* print progress information */

    u64 threshold;                /* any integer less than this is a DP */
    u64 dp_max_it;                /* how many iterations to find a DP */
    u64 points_per_version;       /* #DP per version of the function */
    u64 nslots;                   /* entries in the dict */

    static double optimal_theta(double log2_w, int n)
    {
        return std::log2(2.25) + 0.5 * (log2_w - n);
    }

    void finalize(int n, u64 _nslots)
    {
        if (nbytes_memory == 0)
            errx(1, "the size of the dictionnary must be specified");

        nslots = _nslots;
        double log2_w = std::log2(nslots * n_nodes);
        double theta = optimal_theta(log2_w, n);
        /* auto-choose the difficulty if not set */
        if (difficulty < 0) {
            difficulty = -theta;
            if (difficulty < 0)
                difficulty = 0;
            if (verbose)
                printf("AUTO-TUNING: setting difficulty %.2f\n", difficulty);
        } else {
            if (verbose)
                printf("NOTICE: using difficulty=%.2f vs ``optimal''=%.2f\n", difficulty, theta);
        }
        threshold = std::pow(2., n - difficulty);
        dp_max_it = 20 * std::pow(2., difficulty);
        points_per_version = beta * nslots;

        /* display warnings if problematic choices were made */
        if (verbose && difficulty <= 0) {
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
            printf("---> zero difficulty (use the naive technique!)\n");
            printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");            
        }
        if (verbose && n - difficulty < log2_w) {
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

/* try to iterate for 1s. Return #it/s */
template<typename ConcreteProblem>
static double sequential_benchmark(const ConcreteProblem& Pb)
{
    u64 i = 42;
    u64 threshold = (1ull << (Pb.n - 1));
    double start = wtime();
    u64 total = 0;
    u64 k = 1000;
    u64 x = 0;
    double delta = 0.;
    
    while (delta < 1.) {
        for (u64 j = 0; j < k; j++) {
            u64 y = Pb.mixf(i, x);
            if (is_distinguished_point(y, threshold))
                x = j;
            else
                x = y;
        }
        total += k;
        delta = wtime() - start;
    }
    return total / delta;
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


template<class ConcreteProblem>
std::optional<std::tuple<u64,u64,u64>> process_distinguished_point(
    const ConcreteProblem &Pb, Counters &ctr, Dict<std::pair<u64, u64>> &dict, 
    u64 i, u64 start0, u64 end, u64 len0)
{
    auto probe = dict.pop_insert(end, std::pair(start0, len0));
    if (not probe) {
        ctr.probe_failure();
        return std::nullopt;
    }

    auto [start1, len1] = *probe;
    auto collision = walk(Pb, i, start0, len0, start1, len1);  
    if (not collision) {
        ctr.walk_failure();
        return std::nullopt;         /* robin-hood, or dict false positive */
    }

    auto [x0, x1, y] = *collision;
    if (x0 == x1) {
        ctr.collision_failure();
        return std::nullopt;    /* duh */
    }

    ctr.found_collision();
    if (Pb.mix_good_pair(i, x0, x1)) {
        printf("\nFound golden collision!\n");
        return std::make_optional(std::tuple(i, x0, x1));
    }
    return std::nullopt;
}


}
#endif