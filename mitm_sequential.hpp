#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <iostream>
#include <assert.h>




/// Concepts to enforce on the type. This requires C++20 support. If you use intel compiler,
/// it has to be a version above or equal ICX 2023.1 (oneAPI 2023.1)
/* Ensures that a type B is derived from type B */
//template <typename T>
//concept Equalable = requires (T& a, T& b){
//    { a == b } -> std::convertible_to<bool>;
//};
//
//template <typename T>
//concept Assignable = requires (T a, T b){
//    { a = b } ;
//};


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
// requires Equalable<repr> && Assignable<repr> /* repr must support equality-testing and assignment */
template<typename repr>
class AbstractDomain {
public:
    using t = repr;            /* t is the machine representation of elements of the domain */
    template<typename PRNG>
    void randomize(t &x, PRNG &p) const;   /* set x to a random value */

    void next(t &p); /* start with an element p, go to next element and write that to p */

    int length;                                    /* size in bytes of the serialization */
    size_t n_elements;                             /* How many elements are there in this domain */
    void serialize(const t &x, void *out) const;   /* write this to out */
    void unserialize(t &x, void *in) const;        /* read this from in */

    /* Either we make the output length as an argument or use two function */
    auto checksum(const t &x) const -> uint64_t ;  /* return 1 bit from this */
    auto hash(const t &x) const -> uint64_t ;      /* return more bits from this */

};

/* Ensures that a type Instance is an instance of a AbstractProblem class  */
/// Source: https://stackoverflow.com/a/71921982
//template <typename Instance>
//concept InstanceOfAbstractDomain = requires(Instance c) {
//    // IILE, that only binds to A<...> specialisations
//    // Including classes derived from them
//    []<typename X>( AbstractDomain<X>& ){}(c);
//};

/*
 * Provides the description of the type A and a function f : A -> A.
 *
 * The problem instance can contain extra data.
 * E.g. In the attack on double-encryption, where the goal is to find 
 *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
 */
// requires Equalable<repr> && Assignable<repr> /* repr must support equality-testing and assignment */
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
 * The main engine of the code
 */
template<typename Pb>
auto collision(Pb &pb) -> std::pair<typename Pb::A::t, typename Pb::A::t>
{
    // planning:
    // 1- create a table with all with images of all elements in A, i.e. f(A)
    // 2- sort this table
    // 3- create another table to store collisions
    // 4- walk again through all the images of A under g, i.e. g(A)
    // 5- each time test for collision.

    using t = typename Pb::A::t;

    t v1, v2; // dummy variables



    // create array for type that contains all outputs of f with inputs
    // sort it
    // walk along

    // enforce that Pb is a subclass of AbstractProblem
    static_assert(std::is_base_of<AbstractProblem<typename Pb::A>, Pb>::value, 
            "Pb not derived from AbstractProblem");
    return std::pair(v1, v2);
}





#endif