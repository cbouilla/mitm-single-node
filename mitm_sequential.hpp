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
    static void copy(t& out, t& inp); /* deepcopy inp to out */
    inline static auto extract_1_bit(t& inp) -> int;
    inline static auto extract_k_bits(t& inp, int k) -> uint64_t;

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

    static void send_C_to_A(A_t& out_A, C_t& inp_C);
    static void send_C_to_B(B_t& out_B, C_t& inp_C);
};



/// Since we check the dist point at the beginning we don't need to use the length explicitly.
/// Thus, no need write two version of the function below (one starts with f, the other starts with g).
///  Instead, switch the order of the arguments
template<typename Problem, typename C_t>
auto generate_dist_point(const int theta, /* how many bits should be zero */
                         C_t& inp_A, /* WARNING: this function will edit this argument */
                         C_t& inp_B,
                         C_t& out_C,
                         uint8_t* out_C_serialized)
{

    // void (*f)(A_t&, C_t& ), /* why don't we pass the problem instead */
    // void (*g)(B_t&, C_t& ),
    // void (*send_C_to_A)(A_t&, C_t),
    // void (*send_C_to_B)(B_t&, C_t),
    // void (serialize)(const C_t&, uint8_t*),

    /// inp = 1bit(f/g) || dist point ||
    // What is the user passed c-style array as A_t. Abandon this idea.
    // static A_t inp_A = inp_A_orig; /* make sure this is a deepcopy */
    //static B_t inp_B; /* since we may iterate using g : B -> C  */
    bool found_distinguished = false;
    int f_or_g; // 1 if the next function is f,


    // Do 1 round without entering the loop since we know which function to use
    f(inp_A, out_C);
    // serialize(out_C, out_C_serialized);
    /* decide what is the next function based on the output */
    f_or_g = extract_1_bit(out_C);// next_f_or_g(out_C_serialized);
    // remove the bit used to decide which function


    out_C_serialized[0] = out_C_serialized[0]>>1;
    found_distinguished = ;

    /* potentially infinite loop, todo limit  the number of iteration as a function of theta */
    while (not found_distinguished){

        if (f_or_g == 1){ // i.e. next use f to iterate
            // summary:  C -> A -f-> C
            /* convert output to A input */
            send_C_to_A(inp_A, out_C);
            f(inp_A, out_C);
        } else { // use g in the sequence
            // summary:  C -> B -g-> C
            send_C_to_B(inp_B, out_C);
            f(inp_A, out_C);
        }

        serialize(out_C, out_C_serialized);
        /* decide what is the next function based on the output */
        f_or_g = next_f_or_g(out_C_serialized);
        // remove the bit used to decide which function
        out_C_serialized[0] = out_C_serialized[0]>>1;
        found_distinguished = is_distinguished(out_C_serialized, theta);
    }

    /* For later use when we are going to bound the loop */
    return true;
}


template<typename A_t, typename B_t, typename C_t>
auto generate_seq(void (*f)(A_t&, C_t& ),
                  void (*g)(B_t&, C_t& ),
                  void (*send_C_to_A)(A_t&, C_t),
                  void (*send_C_to_B)(B_t&, C_t),
                  void (serialize)(const C_t&, uint8_t*),
                  C_t& inp1,
                  C_t& inp2,
                  const int theta)
                  -> int

{
    /// Assume that |C| = |A| = |B| for now, this assumption will change !
    std::vector<C_t> inp1_array(3*theta); /* inp 1 output chain */
    std::vector<C_t> inp2_array(3*theta); /* inp 2 output chain */
    uint8_t inp1_serialized[C::length];


    int found_dist_inp1 = 0;
    size_t chain_length_inp1 = 0;
    int found_dist_inp2 = 0;
    size_t chain_length_inp2 = 0;
    int f_or_g;

    for (; chain_length_inp1 < 3*theta; ++chain_length_inp1 ){
        /* iterate */
        f_or_g = next_f_or_g(inp1)

        if (found_dist_inp1)
            break; // exit the loop
    }


}

template<typename A_t, typename B_t, typename C_t>
auto walk(void (*f)(A_t&, C_t& ),
          void (*g)(B_t&, C_t& ),
          void (*send_C_to_A)(A_t&, C_t),
          void (*send_C_to_B)(B_t&, C_t),
          void (serialize)(const C_t&, uint8_t*),
          C_t& inp1,
          C_t& inp2,
          const int theta)
          -> std::pair<A_t, B_t>
{
    /// Given two inputs of different types walk along the two sequences till you find a common point.
    /// This is done by first going through the first sequence till a dist point is found.
    /// Repeat the same with second sequence.
    /// Walk backward in the two computed sequences until you find a uncommon point, stop there.
    /// add a drawing to illustrate this.
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
    using t_pair = typename std::pair<A_t, C_t>;

    // enforce that Pb is an subclass of AbstractProblem
    // todo: rewrite this
    // static_assert(std::is_base_of<AbstractProblem<Domain_A>, Pb>::value,
    //               "Pb not derived from AbstractProblem");


    // -----------------------------------------------------------------------------/
    // --------------------------- INIT --------------------------------
    // DICT
    size_t n_slots = 100;
    Dict<C_t> dict{n_slots}; /* create a dictionary */

    // Pseudo-Random Number Generator



    // -----------------------------------------------------------------------------/
    // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
    int theta = 1; // difficulty;
    A_t a{}; /* input  */
    A_t a_tmp{};
    B_t b{}; /* input */
    B_t b_tmp{};
    C_t c_1{}; /* output */
    C_t c_2{};
    uint8_t c_serial_1[Domain_C::length];
    uint8_t c_serial_2[Domain_C::length];
    // sample random inputs
    C_t c{}; /* main input   */
    /* fill the input */
    Domain_C::randomize(c, pr); // todo add PRG
    uint8_t c_serial[Domain_C::length];
    Domain_C::serialize(c, c_serial);
    int f_or_g = c_serial[0]&1; // 1 if we use f, use g otherwise.


    // loop until enough collisions is collected.
    size_t n_collisions = 0;
    size_t n_needed_collisions = 1;

    while (n_collisions < n_needed_collisions){
        if (f_or_g){
            Domain_A::copy(a_tmp, a); /* copy(out, inp) */
            generate_dist_point(&Pb::f,
                                &Pb::g,
                                &Pb::send_C_to_A,
                                &Pb::send_C_to_B,
                                &Domain_C::serialize,
                                theta,
                                a_tmp,
                                b_tmp,
                                c_1,
                                c_serial_1);

        } else {
            Domain_B::copy(b_tmp, b); /* copy(out, inp) */
            generate_dist_point( &Pb::g, /* reversed order */
                                 &Pb::f,
                                 &Pb::send_C_to_B,
                                 &Pb::send_C_to_A,
                                 &Domain_C::serialize,
                                 theta,
                                 b_tmp,
                                 a_tmp,
                                 c_2, /* why don't we pass b as well? */
                                 c_serial_2);
        }


    }

    // send them to dict
    /* store all pairs of (input, output) in an array. The inverse order to sort in the first element */
    std::vector< t_pair > all_images(dom_A.n_elements);




    return std::pair(a, b); // dummy value
}
#endif