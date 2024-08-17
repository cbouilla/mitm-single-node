#ifndef MITM
#define MITM

#include <unordered_map>

#include "tools.hpp"
#include "problem.hpp"

/*
 * Sequential very naive MITM.
 */

namespace mitm {

template <class AbstractProblem>
optional<pair<u64, u64>> naive_collision_search(AbstractProblem &Pb)
{
    static_assert(std::is_base_of<AbstractCollisionProblem, AbstractProblem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");

    u64 N = 1ull << Pb.n;
    std::unordered_multimap<u64, u64> A;
    A.reserve(N);
    for (u64 x = 0; x < N; x++) {
        u64 z = Pb.f(x);
        A.emplace(z, x);
    }

    for (u64 y = 0; y < N; y++) {
        u64 z = Pb.f(y);
        auto range = A.equal_range(z);
        for (auto it = range.first; it != range.second; ++it) {
            u64 x = it->second;
            if (x != y && Pb.is_good_pair(x, y))
                return std::make_optional(pair(x, y));
        }
    }
    return std::nullopt;
}


template <class AbstractProblem>
vector<pair<u64, u64>> naive_claw_search(AbstractProblem &Pb)
{
    static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
        "problem not derived from mitm::AbstractClawProblem");
  
    double start = wtime();
    u64 N = 1ull << Pb.n;
    std::unordered_multimap<u64, u64> A;
    A.reserve(N);
    for (u64 x = 0; x < N; x++) {
        u64 z = Pb.f(x);
        A.emplace(z, x);
    }

    double mid = wtime();
    printf("Fill: %.1fs\n", mid - start);
    vector<pair<u64, u64>> result;
    for (u64 y = 0; y < N; y++) {
        u64 z = Pb.g(y);
        auto range = A.equal_range(z);
        for (auto it = range.first; it != range.second; ++it) {
            u64 x = it->second;
            if (Pb.is_good_pair(x, y))
                result.push_back(pair(x, y));
        }
    }
    printf("Probe: %.1fs\n", wtime() - mid);
    return result;
}

}

#endif
