#ifndef MITM_TRAIL
#define MITM_TRAIL

#include <cmath>
#include <cassert>
#include <cstdio>

#include "counters.hpp"
#include "parameters.hpp"
#include "problem.hpp"

/*
 * Walking trails: iterate the mixing function to a distinguished point, and turn a
 * dictionary hit into the collision that caused it.
 */

namespace mitm {

inline bool is_distinguished_point(u64 x, u64 threshold)
{
    return x <= threshold;
}

/*
 * Given two inputs that (maybe) lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * This function assumes that the provided lengths (= distance to the distinguished point)
 * are correct, but it does not assume that the two trails end at the same DP.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk(ProblemWrapper& wrapper, Counters &ctr, const Parameters &params,
    u64 i, u64 x0, u64 len0, u64 x1, u64 len1__)
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
    assert(not is_distinguished_point(x1, params.threshold));
    assert(not is_distinguished_point(x0, params.threshold));

    /* move the longest sequence until the remaining number of steps is equal */
    /* to the shortest sequence. */
    u64 len1 = len1__;
    for (; len0 > len1; len0--)
        x0 = wrapper.mixf(i, x0);
    for (; len0 < len1; len1--)
        x1 = wrapper.mixf(i, x1);

    if (x0 == x1) { /* robin-hood */
        ctr.bad_walk_robinhood += 1;
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
            return optional(tuple(x0, x1, len1__));
        }
        x0 = y0;
        x1 = y1;
    }

    if (x0 != x1)    /* false positive from the dictionnary */
        ctr.bad_walk_noncolliding += 1;
    return nullopt;
}

/*
 * Given two inputs that (maybe) lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * This function assumes that len0 (= distance to the distinguished point)
 * is correct, but it does not assume that the two trails end at the same DP.
 * `end` is [end of trail] / params.n_inserters
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> walk_nolen1(ProblemWrapper& wrapper, Counters &ctr, const Parameters &params,
    u64 i, u64 x0, u64 len0, u64 end0, u64 x1)
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

    /* the distance from x1 to a distinguished point is unknown.
     * We need to walk the trail again, but we save all intermediate points.
     */
    u64 maxit = params.dp_max_it;
    u64 trail1[maxit];
    trail1[0] = x1;
    u64 len1 = 0;
    assert(not is_distinguished_point(x1, params.threshold));
    assert(not is_distinguished_point(x0, params.threshold));
    for (;;) {
        len1 += 1;
        x1 = wrapper.mixf(i, x1);
        trail1[len1] = x1;
        if (is_distinguished_point(x1, params.threshold))
            break;
    }

    if (x1 / params.n_inserters != end0) {
        ctr.bad_walk_noncolliding += 1;
        return nullopt;
    }

    /* move the longest sequence until the remaining number of steps is equal */
    /* to the shortest sequence. */
    for (; len0 > len1; len0--)
        x0 = wrapper.mixf(i, x0);

    /* at this stage, len0 <= len1 */
    x1 = trail1[len1 - len0];
    if (x0 == x1) { /* robin-hood */
        ctr.bad_walk_robinhood += 1;
        return nullopt;
    }

    /* now both sequences needs exactly `len0` steps to reach the common distinguished point */
    for (u64 j = len1 - len0;; j++) {
        /* walk them together */
        u64 y0 = wrapper.mixf(i, x0);
        u64 y1 = trail1[j+1];
        /* do the outputs collide? If yes, return true and exit. */
        if (y0 == y1) {
            /* careful: x0 & x1 contain inputs before mixing */
            return optional(tuple(x0, x1, len1));
        }
        x0 = y0;
        x1 = y1;
    }
}

/*
 * Walker side of the engine: given a dictionary hit, walk both trails to locate the
 * collision and test whether it is the golden pair.  Returns (i, x0, x1) if it is.
 *
 * This is the expensive half of processing a distinguished point, which is exactly
 * why it does not run on the inserter thread that found the hit.
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> resolve_collision(ProblemWrapper &wrapper, Counters &ctr, const Parameters &params,
                                               u64 i, u64 root_seed, u64 seed0, u64 end, u64 len0,
                                               u64 seed1, u64 len1_maybe)
{
    u64 start0 = (root_seed + params.multiplier * seed0) & wrapper.out_mask;
    u64 start1 = (root_seed + params.multiplier * seed1) & wrapper.out_mask;

    optional<tuple<u64,u64,u64>> collision;
    if (len1_maybe == 0)
        collision = walk_nolen1(wrapper, ctr, params, i, start0, len0, end, start1);
    else
        collision = walk(wrapper, ctr, params, i, start0, len0, start1, len1_maybe);

    if (not collision)
        return nullopt;         /* robin-hood, or dict false positive */

    auto [x0, x1, len1] = *collision;
    assert(len1_maybe == 0 || len1_maybe == len1);
    if (x0 == x1) {
        ctr.bad_collision += 1;
        return nullopt;    /* duh */
    }

    u64 y0 = wrapper.mix(i, x0);
    u64 y1 = wrapper.mix(i, x1);
    assert(wrapper.mixf(i, x0) == wrapper.mixf(i, x1));
    ctr.found_collision(std::min(y0, y1), len0, std::max(y0, y1), len1);

    if (wrapper.mix_good_pair(i, x0, x1)) {
        printf("\nFound golden collision! i=%" PRIx64 " root_seed=%" PRIx64 " seed0=%" PRIx64 ". Dict --> seed1=%" PRIx64 "\n",
            i, root_seed, seed0, seed1);
        return optional(tuple(i, x0, x1));
    }
    return nullopt;
}

}
#endif
