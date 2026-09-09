#ifndef MITM_SHA2_PROBLEM
#define MITM_SHA2_PROBLEM

#include <cassert>

#include "problem.hpp"

/* defined in sha256.c */
extern "C" {
    void sha256_process(u32 state[8], const u8 data[], u32 length);
}

namespace mitm {

static const u32 sha256_IV[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

/*
 * One compression-function call on a 64-byte block whose first 8 bytes are `x` and
 * whose remaining bytes are all `filler`, truncated to n bits.  This is the common
 * half of both problems below; `filler` is what makes f and g different.
 */
static u64 sha2_truncated(u64 x, u64 filler, u64 mask)
{
    u32 state[8];
    for (int i = 0; i < 8; i++)
        state[i] = sha256_IV[i];
    u64 msg[8];
    for (int i = 0; i < 8; i++)
        msg[i] = filler;
    msg[0] = x;
    sha256_process(state, (const u8*) msg, 64);
    return (state[0] ^ ((u64) state[1] << 32)) & mask;
}


/*
 * Find x != y with f(x) == f(y), where f is SHA256's compression function truncated
 * to n bits.  A collision is planted by making f map golden_y onto f(golden_x), so
 * that the demo terminates on a problem size we can actually run.
 */
class SHA2CollisionProblem : public AbstractCollisionProblem {
public:
    int n, m;
    u64 mask;
    u64 golden_x, golden_y;         /* cheating */

    SHA2CollisionProblem(int n, PRNG &prng) : n(n), m(n)
    {
        mask = make_mask(n);
        golden_x = prng.rand() & mask;
        golden_y = prng.rand() & mask;
        while (golden_y == golden_x)
            golden_y = prng.rand() & mask;
        assert(f(golden_x) == f(golden_y));
    }

    u64 f(u64 x) const
    {
        assert((x & mask) == x);
        if (x == golden_y)
            x = golden_x;
        return sha2_truncated(x, 0, mask);
    }

    /* symmetric, as the contract requires: an engine that finds the pair in either order retires it */
    bool is_good_pair(u64 x0, u64 x1) const
    {
        return ((x0 == golden_x) && (x1 == golden_y)) || ((x0 == golden_y) && (x1 == golden_x));
    }
};


/*
 * Find x, y with f(x) == g(y), where f and g are SHA256's compression function on two
 * different blocks, truncated to n bits.  Again the golden claw is planted, this time
 * by xoring a constant into g.
 */
class SHA2ClawProblem : public AbstractClawProblem {
public:
    int n, m;
    u64 mask;
    u64 g_shift, golden_x, golden_y;    /* cheating */

    SHA2ClawProblem(int n, PRNG &prng) : n(n), m(n)
    {
        mask = make_mask(n);
        golden_x = prng.rand() & mask;
        golden_y = prng.rand() & mask;
        g_shift = 0;
        g_shift = g(golden_y) ^ f(golden_x);
        assert(f(golden_x) == g(golden_y));
    }

    u64 f(u64 x) const
    {
        assert((x & mask) == x);
        return sha2_truncated(x, 0, mask);
    }

    u64 g(u64 x) const
    {
        assert((x & mask) == x);
        return sha2_truncated(x, 0xffffffff, mask) ^ g_shift;
    }

    bool is_good_pair(u64 x, u64 y) const
    {
        return (x == golden_x) && (y == golden_y);
    }
};

}
#endif
