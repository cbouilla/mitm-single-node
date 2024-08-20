#ifndef MITM
#define MITM

#include <optional>
#include <unordered_map>
#include <cassert>

#include "common.hpp"
#include "problem.hpp"
#include "engine_common.hpp"

namespace mitm {

template <typename AbstractProblem>
class ConcreteCollisionProblem {
public:
    const AbstractProblem &pb;
    const int n;
    const u64 mask;
    u64 n_eval;      // #evaluations of (mix)f.  This does not count the invocations of f() by pb.good_pair().

    ConcreteCollisionProblem(const AbstractProblem &pb) : pb(pb), n(pb.n), mask((1ull << pb.n) - 1)
    {
        static_assert(std::is_base_of<AbstractCollisionProblem, AbstractProblem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");
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

    auto [i, x, y] = _Engine::run(wrapper, params, prng);
    u64 a = wrapper.mix(i, x);
    u64 b = wrapper.mix(i, y);
    assert(a != b);
    assert(Pb.f(a) == Pb.f(b));
    assert(Pb.is_good_pair(a, b));
    return pair(a, b);
}

/****************************************************************************************/

template <class AbstractProblem>
class ClawWrapper {
public:
    const AbstractProblem& pb;
    const int n;
    const u64 mask;
    u64 n_eval;        // #evaluations of (mix)f.  This does not count the invocations of f() by pb.good_pair().

    ClawWrapper(const AbstractProblem& pb) : pb(pb), n(pb.n), mask((1ull << pb.n) - 1)
    {
        static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
            "problem not derived from mitm::AbstractClawProblem");
    }

    /* pick either f() or g() */
    bool choose(u64 i, u64 x) const
    {
        return ((x * (i | 1))  >> (n - 1)) & 1;
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

template <typename _Engine, class Parameters, typename Problem>
pair<u64, u64> claw_search(const Problem& Pb, Parameters &params, PRNG &prng)
{
    static_assert(std::is_base_of<Engine, _Engine>::value,
            "engine not derived from mitm::Engine");

    ClawWrapper<Problem> wrapper(Pb);

    auto [i, a, b] = _Engine::run(wrapper, params, prng);
    auto [u, v] = wrapper.swap(i, a, b);
    u64 x0 = wrapper.mix(i, u);
    u64 x1 = wrapper.mix(i, v);
    assert(Pb.is_good_pair(x0, x1));
    return pair(x0, x1);
}


}

#endif
