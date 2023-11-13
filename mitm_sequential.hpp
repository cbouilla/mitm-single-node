#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "include/dict.hpp"

/******************************************************************************/
/* Document for standard implementation                                       */
/******************************************************************************/

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
    static int length;  /* nbytes needed to encode an element */                                   /* size in bytes of the serialization */
    static size_t n_elements; /* how many elements in the domain */
    using t = repr;            /* t is the machine representation of elements of the domain */

    template<class PRNG>
    static void randomize(t &x, PRNG &p);           /* set x to a random value */

    static void randomize(t &x);  /* set x to a random value */



    /* get the next element after x. What matters is getting a different element each time, not the order. */
    auto next(t& x) -> t;
    static void serialize(const t &x, const uint8_t *out);   /* write this to out */
    static void unserialize(t &x, const uint8_t *in);        /* read this from in */
    static void copy(t& out, const t& inp); /* deepcopy inp to out */
    inline static auto extract_1_bit(const t& inp) -> int;
    inline static auto extract_k_bits(const t& inp, int k) -> uint64_t;

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
    /* these lines have to be retyped again */
    using C =  Domain_C;
    using C_t = typename C::t;

    using A = Domain_A;
    using A_t = typename A::t;


    using B = Domain_B;
    using B_t = typename B::t;

    static void f(const A_t &x, C_t &y);                /* y <--- f(x) */
    static void g(const B_t &x, C_t &y);                /* y <--- g(x) */


    AbstractProblem() {
        // enforce that A is a subclass of AbstractDomain
        static_assert(std::is_base_of<AbstractDomain<typename A::t>, A>::value,
                      "A not derived from AbstractDomain");
    }

    static void send_C_to_A(A_t& out_A, C_t& inp_C);
    static void send_C_to_B(B_t& out_B, C_t& inp_C);
};

/******************************************************************************/
/* a user does not need to look at the code below                             */
/******************************************************************************/
template <typename P>
struct Iterate_F {
    ///  A wrapper for calling f that uses
    using A_t = typename P::A::t;
    using C_t = typename P::C::t;
    inline static A_t inp_A; /* A placeholder for the input */
    Iterate_F() {}

    /* Make this struct callable */
    void operator()(C_t& inp_C, C_t& out_C){
        /* convert inp:C_T -> inp:A_t */
        /* enhancement: in the future we can make this function depends on PRNG */
        /* that is defined within this struct. So, we can change the function easily */
        P::send_C_to_A(inp_A, inp_C);
        P::f(inp_A, out_C);
    }
};


template <typename Problem>
struct Iterate_G : Problem{
     /* accessing Problem::g will generate an error*/
    // Problem::g; /* Original iteration function */
    // using Problem::send_C_to_B;
    using B_t = typename Problem::B::t;
    using C_t = typename Problem::C::t;

    // todo should not be local to each thread?
    inline static B_t inp_B{}; /* A placeholder for the input */

    Iterate_G() {}

    /* Make this struct callable */
    void operator()(C_t& inp_C, C_t& out_C){

        /* convert inp:C_T -> inp:A_t */
        /* enhancement: in the future we can make this function depends on PRNG */
        /* that is defined within this struct. So, we can change the function easily */
        Problem::send_C_to_B(inp_B, inp_C);
        Problem::g(inp_B, out_C);
    }
};



/// Since we check the dist point at the beginning we don't need to use the length explicitly.
/// Thus, no need write two version of the function below (one starts with f, the other starts with g).
///  Instead, switch the order of the arguments
template<typename Problem>
auto generate_dist_point(const int64_t theta, /* #bits are zero at the beginning */
                         typename Problem::C::t& inp_C, /* WARNING: this function will edit this argument */
                         typename Problem::C::t& out_C,
                         uint8_t* out_C_serialized,
                         Iterate_F<Problem>& F,
                         Iterate_G<Problem>& G)
                         -> bool
{

    static const uint64_t mask = (1LL<<theta) - 1;

    // using C_t = typename Problem::C::t;
    using C = typename  Problem::C;
    /* F: Domain_C -> Domain_C, instead of C_t -> A_t -f-> C_t */

    bool found_distinguished = false;
    int f_or_g; // 1 if the next function is f,


    /* The decision is based on the input for the next iteration, now it's inp_C then out_C */
    f_or_g = C::extract_1_bit(inp_C);// next_f_or_g(out_C_serialized);

    found_distinguished = (0 == (mask & C::extract_k_bits(inp_C, theta))  );

    /* potentially infinite loop, todo limit  the number of iteration as a function of theta */
    while (not found_distinguished){
        // tood bug switch input and output
        if (f_or_g == 1){ // i.e. next use f to iterate
            // summary:  Domain_C -> A -f-> Domain_C
            /* convert output to A input */
            // send_C_to_A(inp_A, out_C);
           //  f(inp_A, out_C);
           F(inp_C, out_C);
        } else { // use g in the sequence
            // summary:  Domain_C -> Domain_B -g-> Domain_C
            // send_C_to_B(inp_B, out_C);
            // g(inp_B, out_C);
            G(inp_C, out_C);
        }

        //
        /* decide what is the next function based on the output */
        f_or_g = C::extract_1_bit(out_C);
    }

    /* save the output in a serial form  */
    C::serialize(out_C, out_C_serialized);
    /* For later use when we are going to bound the loop */
    return true;
}


template<typename Problem>
auto fill_sequence(typename Problem::C::t& inp_C,
                   std::vector<typename Problem::C::t>& inp_array,
                   const int theta,
                   Iterate_F<Problem>& F,
                   Iterate_G<Problem>& G)
                   -> size_t
{
    /// Given an input C_t inp, fill the inp_array with all output of f/g (inp)
    /// until a distinguished point is found. Return the number of steps needed
    /// to arrive at a distinguished point.
    using C = typename Problem::C;
    using C_t = typename Problem::C::t;


    static C_t out_C{};
    const uint64_t mask = (theta<<1LL) - 1;
    size_t chain_length = 0;
    int found_dist = 0;
    int f_or_g = C::extract_1_bit(inp_C);

    for (; chain_length < (3 * (theta << 1LL)) ; chain_length += 2 ){
        /* do 2 rounds at once to avoid input swapping with output */
        /* iterate: 1 out of 2  */
        if (f_or_g){
            F(inp_C, out_C);
        } else {
            G(inp_C, out_C);
        }

        f_or_g = C::extract_1_bit(out_C);
        found_dist = (0 == (mask & C::extract_k_bits(out_C, theta)));

        if (found_dist) break; // exit the loop

        /* iterate: 2 out of 2  */
        if (f_or_g){
            F(inp_C, out_C);
        } else {
            G(inp_C, out_C);
        }

        f_or_g = C::extract_1_bit(out_C);
        found_dist = (0 == (mask & C::extract_k_bits(out_C, theta)));

        if (found_dist) {
            ++chain_length; /* one step more */
            break; // exit the loop
        }
    }

    return chain_length;
}




template<typename Problem>
auto walk(typename Problem::C::t& inp1,
          typename Problem::C::t& inp2,
          const uint64_t theta,
          Iterate_F<Problem>& F,
          Iterate_G<Problem>& G)
          -> std::pair<typename Problem::C::t, typename Problem::C::t>
{
    /// Given two inputs that lead to the same distinguished point,
    /// find the earliest collision in the sequence before the distinguished point
    /// add a drawing to illustrate this.

    using C_t = typename Problem::C::t;
    std::vector<C_t> inp2_array(3*theta); /* inp 2 output chain */
    std::vector<C_t> inp1_array(3*theta); /* inp 1 output chain */


    size_t inp1_chain_length = fill_sequence(inp1,
                                             inp1_array,
                                             theta,
                                             F, /* Iteration function */
                                             G); /* Iteration function */

    size_t inp2_chain_length = fill_sequence(inp2,
                                             inp2_array,
                                             theta,
                                             F, /* Iteration function */
                                             G); /* Iteration function */


    /* now walk backward until you find the last point where they share the */
    size_t i = 1;
    for (; i < 3*(1LL<<theta); ++i) {
        if (inp1_array[inp1_chain_length - i] != inp2_array[inp2_chain_length]){
            break;
        }
    }
    size_t idx_inp1 = inp1_chain_length - i;
    size_t idx_inp2 = inp2_chain_length - i;

    /* I don't have a good feeling about the line below */
    return std::pair(inp1_array[idx_inp1], inp2_array[idx_inp2]);
}

template <typename Problem >
auto treat_collision(typename Problem::C::t& inp1,
                     typename Problem::C::t& inp2,
                     const uint64_t theta,
                     std::vector< std::pair<typename Problem::C::t, typename Problem::C::t> >& container,
                     Iterate_F<Problem>& F,
                     Iterate_G<Problem>& G)
                     -> bool {
    /* Convert the input types to A_t and B_t. */
    /*  Either save it to disk then later */
    using A_t = typename Problem::A::t;
    using B_t = typename Problem::B::t;
    using C_t = typename Problem::C::t;

    std::pair<C_t, C_t> pair = walk<Problem>(inp1, /* and inp2 leads to the same dist point */
                                             inp2,
                                             theta, /* #zero bits at the beginning */
                                             F, /* Iteration function */
                                             G); /* Iteration function */

    /* todo here we should add more tests */
    container.push_back(pair);
    return true;
}


template<typename Problem>
auto collision()
-> std::pair<typename Problem::C::t, typename Problem::C::t>
{
    using A_t = typename Problem::A::t;
    using Domain_A = typename  Problem::A;
    //A dom_A = pb.dom_A;

    using B_t = typename Problem::B::t;
    using Domain_B = typename  Problem::B;
    //A dom_B = pb.dom_B;

    using C_t = typename Problem::C::t;
    using Domain_C = typename  Problem::C;




    /* save some boilerplate typing */
    using t_pair = typename std::pair<A_t, C_t>;

    // ------------------------------------------------------------------------/
    // --------------------------------- INIT --------------------------------/
    // DICT
    size_t n_slots = 1LL<<25; /* base this number on the available memory */
    Dict<C_t, uint32_t> dict{n_slots}; /* create a dictionary */

    // Pseudo-Random Number Generator

    // todo use arrays method to compare the two outputs

    // -----------------------------------------------------------------------------/
    // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
    int theta = 1; // difficulty;
    /* inp/out variables are used as input and output to save one 1 copy */
    C_t inp_C{}; /* input output */
    C_t out_C{}; /* output input */
    C_t tmp_C{}; /* placeholder to save popped values from dict */
    uint8_t c_serial[Domain_C::length];


    /* fill the input */
    Domain_C::randomize(inp_C); // todo how to add an optional PRG


    /* Collisions related variables */
    bool found_collision = false;
    size_t n_collisions = 0;
    size_t n_needed_collisions = 1LL<32;
    /* a:A_t -f-> x <-g- b:B_t */
    std::vector< std::pair<A_t, B_t> >  collisions_container;

    /* Iteration Functions */
    Iterate_F<Problem> F{};
    Iterate_G<Problem> G{};

    while (n_collisions < n_needed_collisions){
        /* Get a distinguished point */
        generate_dist_point<Problem>(theta,
                                          inp_C,
                                          out_C,
                                          c_serial,
                                          F,
                                          G);


        /* send the result to dictionary, check if it has a collision  */
        found_collision = dict.pop_insert(inp_C,
                                          out_C, // todo just wrong
                                          tmp_C, // todo just wrong
                                          Domain_C::extract_k_bits);

        if (found_collision) [[unlikely]]{
            // treat collision
            // todo walk function has not been used!
            treat_collision<Problem>(inp_C,
                                          tmp_C,
                                          theta,
                                          collisions_container,
                                          F,
                                          G);
            ++n_collisions;
        }
    }



    return std::pair(out_C, tmp_C); // todo wrong values
}
#endif
