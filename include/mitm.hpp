#ifndef MITM
#define MITM

#include <cassert>

#include "tools.hpp"
#include "problem.hpp"
#include "parameters.hpp"
#include "engine.hpp"

/*
 * Umbrella header.  A driver includes this one, plus its own problem definition.
 */

namespace mitm {

/*
 * A collision problem as one random function: x -> f(i ^ x).  Needs |Range| >= |Domain|.
 * INCOMPLETE (CLAUDE.md, Gotchas): for some seeds the search plateaus and never retires the golden
 * pair.  Suspected: mix() is a plain xor, so every version i has the same collision set, and
 * mix_good_pair() does not normalise the pair's order the way swapmix() does.
 */
template <class Problem>
class CollisionWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low n bits */
    const u64 out_mask;             /* the low m bits */
    static constexpr int vlen = Problem::vlen;


    CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");
        assert(m <= 64);

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    /* randomization by a family of permutations of {0, 1}^n */
    u64 mix(u64 i, u64 x) const   /* return σ_i(x) */
    {
        return i ^ x;
    }

    /* evaluates f o σ_i(x) */
    u64 mixf(u64 i, u64 x) const
    {
        return pb.f(mix(i, x));
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            y[j] = mix(i, x[j]);
        /* a scalar problem provides no vf(): call f() straight away */
        if constexpr (vlen == 1)
            r[0] = pb.f(y[0]);
        else
            pb.vf(y, r);
    }

    bool mix_good_pair(u64 i, u64 x0, u64 x1) const
    {
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }
};


/****************************************************************************************/

template <class Problem>
class EqualSizeClawWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low n bits */
    const u64 out_mask;             /* the low m bits */
    const u64 choice_mask;          /* the bit of the mixed value that picks f or g */
    static constexpr int vlen = Problem::vlen;

    EqualSizeClawWrapper(const Problem& pb)
        : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)), choice_mask(1ull << (pb.m - 1))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(pb.n == pb.m);

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    /* pick either f() or g() */
    bool choose(u64 i, u64 x) const
    {
        return (x * (i | 1)) & choice_mask;
    }

    u64 mix(u64 i, u64 x) const
    {
        return i ^ x;
    }

    u64 mixf(u64 i, u64 x) const
    {
        u64 y = mix(i, x);
        if (choose(i, x))
            return pb.f(y);
        else
            return pb.g(y);
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choices[vlen];
        for (int j = 0; j < vlen; j++) {
            y[j] = mix(i, x[j]);
            choices[j] = choose(i, x[j]);
        }
        /* a scalar problem provides no vfg(): call f() / g() straight away */
        if constexpr (vlen == 1)
            r[0] = choices[0] ? pb.f(y[0]) : pb.g(y[0]);
        else
            pb.vfg(y, choices, r);
    }

    pair<u64, u64> swapmix(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(mix(i, x0), mix(i, x1));
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) const
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swapmix(i, a, b);
        return pb.is_good_pair(x0, x1);
    }
};

template <class Problem>
class LargerRangeClawWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits: pb.n + 1, the choice bit included */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low pb.n bits */
    const u64 out_mask;             /* the low m bits */
    static constexpr int vlen = Problem::vlen;
    u64 choice_mask;                /* the bit of the mixed value that picks f or g */

    LargerRangeClawWrapper(const Problem& pb)
        : pb(pb), n(pb.n + 1), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(n <= m);
        choice_mask = 1ull << n;

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    inline u64 full_mix(u64 i, u64 x) const
    {
        u64 y = x * (i | 1);
        return (y ^ (y >> n));
    }

    /* pick either f() or g() */
    bool choose(u64 i, u64 x) const
    {
        return full_mix(i, x) & choice_mask;
    }

    u64 mix(u64 i, u64 x) const   // {0, 1}^m  x  {0, 1}^m ---> {0, 1}^n
    {
        return full_mix(i, x) & in_mask;
    }

    u64 mixf(u64 i, u64 x) const  // {0, 1}^m  x  {0, 1}^m ---> {0, 1}^m
    {
        u64 z = full_mix(i, x);
        if (z & choice_mask)
            return pb.f(z & in_mask);
        else
            return pb.g(z & in_mask);
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choice[vlen];
        for (int j = 0; j < vlen; j++) {
            y[j] = mix(i, x[j]);
            choice[j] = choose(i, x[j]);
        }
        /* a scalar problem provides no vfg(): call f() / g() straight away */
        if constexpr (vlen == 1)
            r[0] = choice[0] ? pb.f(y[0]) : pb.g(y[0]);
        else
            pb.vfg(y, choice, r);
    }

    pair<u64, u64> swapmix(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(mix(i, x0), mix(i, x1));
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) const
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swapmix(i, a, b);
        return pb.is_good_pair(x0, x1);
    }
};


/****************************************************************************************/

/* find x0 != x1 with f(x0) == f(x1) and is_good_pair(x0, x1) */
template <class Problem>
optional<pair<u64, u64>> collision_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
    int rank;
    MPI_Comm_rank(opts.mpi_comm, &rank);
    if (opts.verbose && rank == 0)
        printf("Starting collision search with f : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
            pb.n, pb.m, Problem::vlen);

    CollisionWrapper<Problem> wrapper(pb);
    auto collision = run(wrapper, nbytes_memory, opts, prng);
    if (not collision)
        return nullopt;

    auto [i, x, y] = *collision;
    u64 x0 = wrapper.mix(i, x);
    u64 x1 = wrapper.mix(i, y);

    assert(x0 != x1);
    assert(pb.f(x0) == pb.f(x1));
    assert(pb.is_good_pair(x0, x1));
    return optional(pair(x0, x1));
}


/* find x0, x1 with f(x0) == g(x1) and is_good_pair(x0, x1) */
template <class Problem>
optional<pair<u64, u64>> claw_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
    int rank;
    MPI_Comm_rank(opts.mpi_comm, &rank);
    bool verbose = opts.verbose && (rank == 0);
    if (verbose)
        printf("Starting claw search with f, g : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
            pb.n, pb.m, Problem::vlen);

    optional<tuple<u64,u64,u64>> claw;
    u64 x0, x1;

    /* distinct wrapper types, so each branch carries the whole search */
    if (pb.n == pb.m) {
        if (verbose)
            printf("  - using |Domain| == |Range| mode.  Expecting 1.8*n/w rounds.\n");
        EqualSizeClawWrapper<Problem> wrapper(pb);
        claw = run(wrapper, nbytes_memory, opts, prng);
        if (claw) {
            auto [i, a, b] = *claw;
            std::tie(x0, x1) = wrapper.swapmix(i, a, b);
        }
    } else if (pb.n < pb.m) {
        if (verbose)
            printf("  - using |Domain| << |Range| mode.  Expecting 0.9*n/w rounds.\n");
        LargerRangeClawWrapper<Problem> wrapper(pb);
        claw = run(wrapper, nbytes_memory, opts, prng);
        if (claw) {
            auto [i, a, b] = *claw;
            std::tie(x0, x1) = wrapper.swapmix(i, a, b);
        }
    } else {
        errx(1, "Larger domain not yet supported...");
    }

    if (not claw)
        return nullopt;

    assert((x0 & make_mask(pb.n)) == x0);
    assert((x1 & make_mask(pb.n)) == x1);
    assert(pb.f(x0) == pb.g(x1));
    assert(pb.is_good_pair(x0, x1));
    return optional(pair(x0, x1));
}

}
#endif
