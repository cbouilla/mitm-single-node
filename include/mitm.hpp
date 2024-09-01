#ifndef MITM
#define MITM

#include <optional>
#include <unordered_map>
#include <cassert>

#include "common.hpp"
#include "problem.hpp"
#include "engine_common.hpp"

namespace mitm {

/* This works when |Range| >= |Domain| in the input problem */
template <class AbstractProblem>
class ConcreteCollisionProblem {
public:
    const AbstractProblem &pb;
    const int n, m;
    const u64 in_mask, out_mask;
    u64 n_eval;                // #evaluations of (mix)f.  This does not count the invocations of f() by pb.good_pair().


    ConcreteCollisionProblem(const AbstractProblem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractCollisionProblem, AbstractProblem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");
        assert(m <= 64);
        assert(0);   // not reay yet
    }

    /* randomization by a family of permutations of {0, 1}^n */
    u64 mix(u64 i, u64 x) const   /* return σ_i(x) */
    {
        return i ^ x;
    }

    /* evaluates f o σ_i(x) */
    u64 mixf(u64 i, u64 x)
    {
        n_eval += 1;
        return pb.f(mix(i, x));
    }

    bool mix_good_pair(u64 i, u64 x0, u64 x1)
    { 
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }
};


template <typename _Engine, class Parameters, typename AbstractProblem>
pair<u64, u64> collision_search(const AbstractProblem& Pb, Parameters &params, PRNG &prng)
{
    static_assert(std::is_base_of<Engine, _Engine>::value,
            "engine not derived from mitm::Engine");

    ConcreteCollisionProblem wrapper(Pb);

    params.finalize(Pb.n, Pb.m);
    auto [i, x, y] = _Engine::run(wrapper, params, prng);
    u64 a = wrapper.mix(i, x);
    u64 b = wrapper.mix(i, y);
    assert(a != b);
    assert(Pb.f(a) == Pb.f(b));
    assert(Pb.is_good_pair(a, b));
    return pair(a, b);
}

/****************************************************************************************/

// code deduplication could be achieved with the CRTP...

template <class Problem>
class EqualSizeClawWrapper {
public:
    const Problem &pb;
    const int n, m;
    const u64 in_mask, out_mask, choice_mask;
    static constexpr int vlen = Problem::vlen;
    u64 n_eval;                // #evaluations of (mix)f.  This does not count the invocations of f() by pb.good_pair().

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
        u64 x[pb.vlen] __attribute__ ((aligned(sizeof(u64) * pb.vlen))); 
        u64 y[pb.vlen] __attribute__ ((aligned(sizeof(u64) * pb.vlen)));
        for (int i = 0; i < pb.vlen; i++)
            x[i] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < pb.vlen; j++)
            assert(y[j] == mixf(i, x[j]));
        n_eval = 0;
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

    u64 mixf(u64 i, u64 x)
    {
        n_eval += 1;
        u64 y = mix(i, x);
        if (choose(i, x))
            return pb.f(y);
        else
            return pb.g(y);
    }

    void vmixf(u64 i, u64 x[], u64 r[])
    {
        // careful: vlen can be more than one SIMD vector
        constexpr int vlen = Problem::vlen; 
        n_eval += vlen;
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen))); 
        u64 fy[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 gy[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            y[j] = mix(i, x[j]);
        pb.vfg(y, fy, gy);
        for (int j = 0; j < vlen; j++)
            r[j] = choose(i, x[j]) ? fy[j] : gy[j];
    }

    pair<u64, u64> swap(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(x0, x1);    
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) 
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swap(i, a, b);
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }
};

template <class Problem>
class LargerRangeClawWrapper {
public:
    const Problem &pb;
    const int n, m;
    const u64 in_mask, out_mask;
    static constexpr int vlen = Problem::vlen;
    u64 n_eval;                // #evaluations of (mix)f.  This does not count the invocations of f() by pb.good_pair().
    u64 choice_mask;

    LargerRangeClawWrapper(const Problem& pb) : pb(pb), n(pb.n + 1), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)) 
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(n <= m);
        choice_mask = 1ull << n;
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

    u64 mixf(u64 i, u64 x)        // {0, 1}^m  x  {0, 1}^m ---> {0, 1}^m
    {
        n_eval += 1;
        u64 z = full_mix(i, x);
        if (z & choice_mask)
            return pb.f(z & in_mask);
        else
            return pb.g(z & in_mask);
    }
    
    void vmixf(u64 i, u64 x[], u64 r[])
    {
        // careful: vlen can be more than one SIMD vector
        constexpr int vlen = Problem::vlen; 
        n_eval += vlen;
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen))); 
        u64 fy[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 gy[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            y[j] = mix(i, x[j]);
        pb.vfg(y, fy, gy);
        for (int j = 0; j < vlen; j++)
            r[j] = choose(i, x[j]) ? fy[j] : gy[j];
    }

    pair<u64, u64> swap(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(x0, x1);    
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) 
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swap(i, a, b);
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }
};


template <class _Engine, class Parameters, class Problem>
pair<u64, u64> claw_search(const Problem& pb, Parameters &params, PRNG &prng)
{
    static_assert(std::is_base_of<Engine, _Engine>::value,
            "engine not derived from mitm::Engine");

    u64 x0, x1;

    if (params.verbose)
        printf("Starting claw search with f : {0,1}^%d --> {0, 1}^%d\n", pb.n, pb.m);

    if (Problem::vlen > 1) {
        if (params.verbose)
            printf("Using vectorized implementation with vectors of size %d\n", pb.vlen);
        // check consistency of the vector function 
        PRNG vprng;
        u64 x[pb.vlen] __attribute__ ((aligned(sizeof(u64) * pb.vlen))); 
        u64 y[pb.vlen] __attribute__ ((aligned(sizeof(u64) * pb.vlen)));
        u64 z[pb.vlen] __attribute__ ((aligned(sizeof(u64) * pb.vlen)));
        u64 mask = make_mask(pb.n);
        for (int i = 0; i < pb.vlen; i++)
            x[i] = vprng.rand() & mask;
        pb.vfg(x, y, z);
        for (int i = 0; i < pb.vlen; i++)
            printf("y[%d] = %" PRIx64 " vs f(x[%d]) = %" PRIx64 "\n", i, y[i], i, pb.f(x[i]));
        for (int i = 0; i < pb.vlen; i++) {
            assert(y[i] == pb.f(x[i]));
            assert(z[i] == pb.g(x[i]));
        }
    }

    if (pb.n == pb.m) {
        if (params.verbose)
            printf("  - using |Domain| == |Range| mode.  Expecting 1.8*n/w rounds.\n");
        EqualSizeClawWrapper<Problem> wrapper(pb);
        params.finalize(wrapper.n, wrapper.m);
        auto [i, a, b] = _Engine::run(wrapper, params, prng);
        auto [u, v] = wrapper.swap(i, a, b);
        x0 = wrapper.mix(i, u);
        x1 = wrapper.mix(i, v);
    } else if (pb.n < pb.m) {
        if (params.verbose)
            printf("  - using |Domain| << |Range| mode.  Expecting 0.9*n/w rounds.\n");
        LargerRangeClawWrapper<Problem> wrapper(pb);
        params.finalize(wrapper.n, wrapper.m);
        auto [i, a, b] = _Engine::run(wrapper, params, prng);
        auto [u, v] = wrapper.swap(i, a, b);
        x0 = wrapper.mix(i, u);
        x1 = wrapper.mix(i, v);
    } else {
        printf("Larger domain not yet supported...\n");
        assert(0);
    }

    /* quality control */
    assert((x0 & ((1 << pb.n) - 1)) == x0);
    assert((x1 & ((1 << pb.n) - 1)) == x1);    
    assert(pb.f(x0) == pb.g(x1));
    assert(pb.is_good_pair(x0, x1));
    return pair(x0, x1);
}


}

#endif
