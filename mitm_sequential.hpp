#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "dict.hpp"


/*
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq
 */
class AbstractPRNG {
public:
    AbstractPRNG(uint64_t seed, uint64_t seq);
    uint64_t rand();
};



/*
 * A "domain" extends some type to provide the extra functions we need.
 * An instance of the domain can contain extra information
 * E.g. for small integers mod p, "repr" could be uint64_t while the domain actually contains p
 * E.g. for points on an elliptic curve, "repr" could be a pair of integers mod p, while the
 *      domain would contain the equation of the curve, etc.
 */
template<class repr>           /* repr must support comparisons, and assignment */
class AbstractDomain {
public:
    using t = repr;            /* t is the machine representation of elements of the domain */
    template<class PRNG>
    static void randomize(t &x, PRNG &p);           /* set x to a random value */

    static int length;                                    /* size in bytes of the serialization */
    static size_t n_elements; /* how many elements in the domain */

    /* get the next element after x. What matters is getting a different element each time, not the order. */
    auto next(t& x) -> t;
    static void serialize(const t &x, void *out);   /* write this to out */
    static void unserialize(t &x, void *in);        /* read this from in */

    static auto hash(const t &x) -> uint64_t ;               /* return some bits from this */
    static auto hash_extra(const t &x) -> uint64_t ;         /* return more bits from this */
};


/*
 * Provides the description of the type A and a function f : A -> A.
 *
 * The problem instance can contain extra data.
 * E.g. In the attack on double-encryption, where the goal is to find
 *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
 */
template<typename Domain_A, typename Domain_B, typename Domain_C>
class AbstractProblem {
public:
    using C = Domain_C; /* output type */
    using C_t = typename C::t;

    using A = Domain_A;
    using A_t = typename A::t;
    static void f(const A_t &x, C_t &y);                /* y <--- f(x) */

    using B = Domain_B;
    using B_t = typename B::t;
    static void g(const B_t &x, C_t &y);                /* y <--- f(x) */


    AbstractProblem() {
        // enforce that Domain_A is a subclass of AbstractDomain
        static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
                      "Domain_A not derived from AbstractDomain");
    }
};


inline
auto is_distinguished(uint8_t* inp, const int theta) -> bool
{
    /// This function is dangerous since it doesn't check that theta is smaller than
    /// than the length of inp! However, we promise that we will only pass theta that
    /// satisifies this condition.
    /// since numbers are stored as small-endian, start with the first byte.
    int result = 0;
    for (int i = 0; i<theta/8; ++i){ // move byte by byte
        /* theta = 8k + i , here we treat the k bytes that must be zero */
        result |= inp[i]; // if we found non-zero value then the result will be non-zero
    }

    /* the remaining bit */
    uint8_t mask =((1<<(8 - (theta % 8))) - 1)<<theta; // 2^(-theta mod 8) - 1, all on
    result |= inp[theta/8] & mask;

    return (result == 0);
}

template<typename A_t, typename B_t, typename C_t>
auto generate_dist_point(void (*f)(A_t&, C_t& ), /* why don't we pass the problem instead */
                         void (*g)(B_t&, C_t& ),
                         void (serialize)(const C_t&, uint8_t*),
                         const int theta, /* how many bits should be zero */
                         A_t& inp_A,
                         C_t& out_C,
                         uint8_t* out_C_serialized)
{
    static B_t out; /* since we may iterate for  */
    bool found_distinguished = false;

    /* First treat the input of A */
    f(inp_A, out_C);
    serialize(out_C, out_C_serialized);
    /* check if it is a distinguished point */

    while (not found_distinguished){

    }

    return;
}

/* yay, function overloading came to rescue! */
template<typename A_t, typename B_t, typename C_t>
auto generate_dist_point(void (*f)(A_t&, C_t& ),
                         void (*g)(B_t&, C_t& ),
                         B_t& inp_B,
                         C_t& out_C)
{
    return;
}





template<typename Pb>
auto collision(const Pb &pb) -> std::pair<typename Pb::A::t, typename Pb::A::t>
{
    using A_t = typename Pb::A::t;
    using Domain_A = typename  Pb::A;
    Domain_A dom_A = pb.dom_A;

    using B_t = typename Pb::B::t;
    using Domain_B = typename  Pb::B;
    Domain_A dom_B = pb.dom_B;

    using C_t = typename Pb::C::t;
    using Domain_C = typename  Pb::C;
    Domain_A dom_C = pb.dom_C;

    /* save some boilerplate typing */
    // todo rename this part :(k
    using t_pair = typename std::pair<A_t, C_t>;

    // enforce that Pb is an subclass of AbstractProblem
    // todo: rewrite this
//    static_assert(std::is_base_of<AbstractProblem<Domain_A>, Pb>::value,
//                  "Pb not derived from AbstractProblem");

    A_t x{}; /* input  */
    C_t y{}; /* output */
    /* store all pairs of (input, output) in an array. The inverse order to sort in the first element */
    std::vector< t_pair > all_images(dom_A.n_elements);

    for (size_t i = 0; i < dom_A.n_elements; ++i) {
        pb.f(x, y);
        /* actual computation */
        all_images[i] = std::pair(x, y); /* let's hope this is a deepcopy */
        dom_A.next(x); /* modifies x */
    }
    /* sort the pairs according to the output */
    std::sort(all_images.begin(),
              all_images.end(),
              [](t_pair x, t_pair y){ return x.second < y.second; }
    );

    t_pair v1 = all_images[0];
    t_pair v2;
    int is_collision_found = 0;
//    /* Check if we have already a collision in the list */
    for (size_t i = 1; i < dom_A.n_elements; ++i){ /* Normally we should not check in the list */
        if (v1.second == all_images[i].second and all_images[i].first != v1.first){
            v2 = all_images[i];
            is_collision_found = 1;
            std::cout << "collision was found in the list!\n";

            return std::pair(v1.first, v2.first); /* break from the function */
            //break;
        }
        v1 = all_images[i]; /* go to the next value */
    }

    /* generate random elements until a collision is found. */
    v1.first = x;
    v1.second = y;
    while(!is_collision_found){
        pb.f(v1.first, v1.second);
        is_collision_found = std::binary_search(all_images.begin(),
                                                all_images.end(),
                                                v1,
                                                [](t_pair p1, t_pair p2){
                                                    return p1.second < p2.second;
                                                }
        );
    }

    // todo this part is rushed
    // find where is the collision then return
    // use std::find to save somelines
    for (size_t i = 0; i < dom_A.n_elements ; ++i) {
        if (all_images[i].second == v1.second and all_images[i].first != v1.first){
            v2 = all_images[i];
            std::cout << "v1.first = " << v1.first << ", v1.second = " << v1.second << "\n" ;
            std::cout << "v2.first = " << v2.first << ", v2.second = " << v2.second << "\n" ;
            break; /* exit the loop */
        }
    }
    return std::pair(v1.first, v2.first);

}
#endif