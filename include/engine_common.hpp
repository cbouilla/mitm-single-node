#ifndef MITM_ENGINE_COMMON
#define MITM_ENGINE_COMMON

#include <cmath>
#include <cassert>
#include <cstdio>

#include "common.hpp"
#include "problem.hpp"
#include "dict.hpp"

namespace mitm {

/******************************************************************************/

class Engine {};

inline bool is_distinguished_point(u64 x, u64 threshold)
{
    return x <= threshold;
}

/*
 * Given an element of the RANGE of f, iterate the function until a distinguished point is found.
 */
template<typename ProblemWrapper>
optional<pair<u64,u64>> generate_dist_point(ProblemWrapper& wrapper, u64 i, const Parameters &params, u64 x)
{
    /* The probability, p, of NOT finding a distinguished point after the loop is
     * Let: theta := 2^-d
     * difficulty, N = k*2^difficulty then,
     * p = (1 - theta)^N =>  let ln(p) <= -k
     */
    for (u64 j = 0; j < params.dp_max_it; j++) {
        u64 y = wrapper.mixf(i, x);
        if (is_distinguished_point(y, params.threshold))
            return optional(pair(y, j + 1));
        x = y;
    }
    return nullopt; /* no distinguished point was found after too many iterations */
}



/* Given two inputs that lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * add a drawing to illustrate this.
 */
template<class ProblemWrapper, class Counters>
optional<tuple<u64,u64,u64>> walk(ProblemWrapper& wrapper, Counters &ctr, u64 i, u64 x0, u64 len0, u64 x1, u64 len1)
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
        x0 = wrapper.mixf(i, x0);
    for (; len0 < len1; len1--)
        x1 = wrapper.mixf(i, x1);
  
    if (x0 == x1) { /* robin-hood */
        ctr.walk_robinhood();
        return nullopt;
    }

    /* now both sequences needs exactly `len` steps to reach the common distinguished point */
    for (u64 j = 0; j < len0; ++j) {
        /* walk them together and check each time if their output are equal     */
        /* return as soon equality is found. */
        u64 y0 = wrapper.mixf(i, x0);
        u64 y1 = wrapper.mixf(i, x1);

        /* First, do the outputs collide? If yes, return true and exit. */
        if (y0 == y1) {
            /* careful: x0 & x1 contain inputs before mixing */
            return std::make_optional(tuple(x0, x1, y0));
        }
        x0 = y0;
        x1 = y1;
    }

    if (x0 != x1)    /* false positive from the dictionnary */
        ctr.walk_noncolliding();
    return nullopt; 
}


template<class ProblemWrapper, class Counters>
optional<tuple<u64,u64,u64>> process_distinguished_point(ProblemWrapper &wrapper, Counters &ctr, PcsDict &dict, 
                                                        u64 i, u64 start0, u64 end, u64 len0)
{
    auto probe = dict.pop_insert(end, start0, len0);
    if (not probe) {
        ctr.probe_failure();
        return nullopt;
    }

    auto [start1, len1] = *probe;
    auto collision = walk(wrapper, ctr, i, start0, len0, start1, len1);  
    if (not collision) 
        return nullopt;         /* robin-hood, or dict false positive */

    auto [x0, x1, y] = *collision;
    if (x0 == x1) {
        ctr.collision_failure();
        return nullopt;    /* duh */
    }

    ctr.found_collision();
    if (wrapper.mix_good_pair(i, x0, x1)) {
        printf("\nFound golden collision!\n");
        return optional(tuple(i, x0, x1));
    }
    return nullopt;
}


}
#endif