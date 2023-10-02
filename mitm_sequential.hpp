#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <iostream>
#include <assert.h>

/* 
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq 
 */
class AbstractPRNG {
public:
    AbstractPRNG(uint64_t seed, uint64_t seq);
    auto rand() -> uint64_t;
};



/*
 * A "domain" extends some type to provide the extra functions we need.
 * An instance of the domain can contain extra information
 * E.g. for small integers mod p, "repr" could be uint64_t while the domain actually contains p
 * E.g. for points on an elliptic curve, "repr" could be a pair of integers mod p, while the
 *      domain would contain the equation of the curve, etc.
 */
template<typename repr>           /* repr must support equality-testing and assignment */
class AbstractDomain {
public:
    using t = repr;            /* t is the machine representation of elements of the domain */
    template<typename PRNG>
    void randomize(t &x, PRNG &p) const;           /* set x to a random value */

    int length;                                    /* size in bytes of the serialization */
    void serialize(const t &x, void *out) const;   /* write this to out */
    void unserialize(t &x, void *in) const;        /* read this from in */

    auto hash(const t &x) const -> uint64_t ;               /* return some bits from this */
    auto hash_extra(const t &x) const -> uint64_t ;         /* return more bits from this */
};


/*
 * Provides the description of the type A and a function f : A -> A.
 *
 * The problem instance can contain extra data.
 * E.g. In the attack on double-encryption, where the goal is to find 
 *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
 */
template<typename Domain>
class AbstractProblem {
public:
    using A = Domain;
    using A_t = typename A::t;
    void f(const A_t &x, A_t &y) const;                /* y <--- f(x) */

    AbstractProblem() {
        // enforce that Domain is a subclass of AbstractDomain
        static_assert(std::is_base_of<AbstractDomain<typename Domain::t>, Domain>::value, 
            "Domain not derived from AbstractDomain");
    }
};

/*
 * Actual collision-finding code (this is Floyd's algorithm)
 */
template<typename Pb>
auto collision(Pb &pb) -> std::pair<typename Pb::A::t, typename Pb::A::t>
{ 
    using t = typename Pb::A::t;

    // enforce that Pb is a subclass of AbstractProblem
    static_assert(std::is_base_of<AbstractProblem<typename Pb::A>, Pb>::value, 
            "Pb not derived from AbstractProblem");
    return ;
}





#endif