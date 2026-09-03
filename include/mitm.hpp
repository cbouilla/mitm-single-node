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
 * Turns a collision problem into one random function: iterate x -> f(i ^ x).
 * Works when |Range| >= |Domain| in the input problem.
 *
 * INCOMPLETE, and known to be so since before the single-engine refactor -- it used
 * to carry a commented-out `assert(0); // not ready yet`.  For some seeds the search
 * plateaus and never retires the golden pair, while the claw wrappers below converge
 * normally.  Verified to behave identically on the deleted sequential engine, so it
 * is the wrapper and not the engine: `mix()` is a plain xor, so every version of the
 * function has the same collision set (unlike the claw wrappers, whose choose() also
 * re-randomizes WHICH of f/g is applied), and mix_good_pair() hands the pair to
 * is_good_pair() in arrival order rather than normalizing it the way swapmix() does.
 */
template <class Problem>
class CollisionWrapper {
public:
    const Problem &pb;
    const int n, m;
    const u64 in_mask, out_mask;
    static constexpr int vlen = Problem::vlen;


    CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");
        assert(m <= 64);

        /* check vmixf against mixf, which also exercises the problem's vf() */
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
        // careful: vlen can be more than one SIMD vector
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

// code deduplication could be achieved with the CRTP...

template <class Problem>
class EqualSizeClawWrapper {
public:
    const Problem &pb;
    const int n, m;
    const u64 in_mask, out_mask, choice_mask;
    static constexpr int vlen = Problem::vlen;

    EqualSizeClawWrapper(const Problem& pb)
        : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)), choice_mask(1ull << (pb.m - 1))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(pb.n == pb.m);

        /* check vmixf */
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
        //u64 z = (i ^ (x * 0xc4ceb9fe1a85ec53ull)) * 0xc6a4a7935bd1e995LLU;
        // return z & in_mask;
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
        // careful: vlen can be more than one SIMD vector
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
    const int n, m;
    const u64 in_mask, out_mask;
    static constexpr int vlen = Problem::vlen;
    u64 choice_mask;

    LargerRangeClawWrapper(const Problem& pb) : pb(pb), n(pb.n + 1), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(n <= m);
        choice_mask = 1ull << n;

        /* check vmixf */
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
        // careful: vlen can be more than one SIMD vector
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

    /* quality control */
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

    /*
     * The two wrappers differ only in how they fold f and g into one function, but
     * they are distinct types, so the branch has to carry the whole search.
     */
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

    /* quality control */
    assert((x0 & make_mask(pb.n)) == x0);
    assert((x1 & make_mask(pb.n)) == x1);
    assert(pb.f(x0) == pb.g(x1));
    assert(pb.is_good_pair(x0, x1));
    return optional(pair(x0, x1));
}

}
#endif
