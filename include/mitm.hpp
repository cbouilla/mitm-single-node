#ifndef MITM
#define MITM

#include <optional>
#include <unordered_map>
#include <cassert>

#include "common.hpp"
#include "problem.hpp"
#include "engine_common.hpp"

namespace mitm {

template <typename AbstractProblem, typename Counters>
class ConcreteCollisionProblem {
public:
    const AbstractProblem &pb;
    Counters &ctr;
    const int n;
    const u64 mask;

    ConcreteCollisionProblem(const AbstractProblem &pb, Counters &ctr) : pb(pb), ctr(ctr), n(pb.n), mask((1ull << pb.n) - 1)
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
    u64 mixf(u64 i, u64 x) const
    {
        ctr.function_evaluation();
        return pb.f(mix(i, x));
    }

  bool mix_good_pair(u64 i, u64 x0, u64 x1) const
  { 
    return pb.is_good_pair(mix(i, x0), mix(i, x1));
  }
};


template <typename Engine, typename AbstractProblem>
pair<u64, u64> collision_search(const AbstractProblem& Pb, Parameters &params, PRNG &prng)
{
    typename Engine::Counters ctr(params.verbose);
    ConcreteCollisionProblem wrapper(Pb, ctr);

    auto [i, x, y] = Engine::run(wrapper, params, prng);
    u64 a = wrapper.mix(i, x);
    u64 b = wrapper.mix(i, y);
    assert(a != b);
    assert(Pb.f(a) == Pb.f(b));
    assert(Pb.is_good_pair(a, b));
    
    ctr.done();
    return pair(a, b);
}

/****************************************************************************************/

template <class AbstractProblem, class Counters>
class ClawWrapper {
public:
    const AbstractProblem& pb;   // AbstractClawProblem
    Counters &ctr;
    const int n;
    const u64 mask;

    ClawWrapper(const AbstractProblem& pb, Counters &ctr) : pb(pb), ctr(ctr), n(pb.n), mask((1ull << pb.n) - 1)
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


    u64 mixf(u64 i, u64 x) const
    {
        ctr.function_evaluation();
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

    bool mix_good_pair(u64 i, u64 a, u64 b) const 
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swap(i, a, b);
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }
};

template <typename Engine, class Parameters, typename Problem>
pair<u64, u64> claw_search(const Problem& Pb, Parameters &params, PRNG &prng)
{
    typename Engine::Counters ctr(params.verbose);
    ClawWrapper wrapper(Pb, ctr);

    auto [i, a, b] = Engine::run(wrapper, params, prng);
    auto [u, v] = wrapper.swap(i, a, b);
    u64 x0 = wrapper.mix(i, u);
    u64 x1 = wrapper.mix(i, v);
    assert(Pb.is_good_pair(x0, x1));

    ctr.done();
    return pair(x0, x1);
}


//=============================================================================+
//----------------------------- NAIVE ENGINES ---------------------------------|

template <typename Problem>
optional<pair<u64, u64>> naive_collision_search(Problem &Pb)
{
  u64 n_items = 1ull << Pb.n;

  std::unordered_multimap<u64, u64> f_images;
  f_images.reserve(n_items);
  for (u64 x = 0; x < n_items; x++) {
    u64 z = Pb.f(x);
    f_images.emplace(z, x);
  }
  assert(f_images.size() == n_items);

  for (u64 y = 0; y < n_items; y++) {
    u64 z = Pb.f(y);
    auto range = f_images.equal_range(z);
    for (auto it = range.first; it != range.second; ++it) {
      u64 x = it->second;
      if (x != y && Pb.is_good_pair(z, x, y))
        return std::make_optional(pair(x, y));
    }
  }
  return std::nullopt;
}


template <typename Problem>
optional<pair<u64, u64>> naive_claw_search(Problem &Pb)
{
  u64 n_items = 1ull << Pb.n;

  std::unordered_multimap<u64, u64> f_images;
  f_images.reserve(n_items);
  for (u64 x = 0; x < n_items; x++) {
    u64 z = Pb.f(x);
    f_images.emplace(z, x);
  }
  assert(f_images.size() == n_items);

  for (u64 y = 0; y < n_items; y++) {
    u64 z = Pb.g(y);
    auto range = f_images.equal_range(z);
    for (auto it = range.first; it != range.second; ++it) {
      u64 x = it->second;
      if (Pb.is_good_pair(z, x, y))
        return std::make_optional(pair(x, y));
    }
  }
  return std::nullopt;
}

}

#endif
