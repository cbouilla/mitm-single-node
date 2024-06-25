#ifndef MITM_ABSTRACT_COLLISION
#define MITM_ABSTRACT_COLLISION
#include "AbstractDomain.hpp"
#include <type_traits>
#include "../include/util/prng.hpp"

namespace mitm {
/*
 * Provides the description of the types A, C, a function f : A -> C.
 * and an optional predicate P.
 *
 * The problem instance can contain extra data, e.g. if the goal consists
 * in finding H(prefix || x) == H(prefix || y) with x != y, then the Problem
 * could contain prefix.
 */
template<typename I, typename A, typename Domain_C>
class AbstractCollisionProblem {
public:
  /* these lines have to be retyped again to use the above 3 types directly */
  using I_t = I; 
  using A_t = A;
  using C_t = typename Domain_C::t;
  using Dom_C = Domain_C;
  Dom_C C;
  /* Used in collecting stats */
  size_t nbits_A;
  size_t nbits_C;

  /* Initialization, in principle  we don't need to initialize any variable unless a user
   * wishes to do so.
   */
  AbstractCollisionProblem() {
    // enforce that C is a subclass of AbstractDomain
    // In fact, we don't need to know **anything** about A
    /* At the end, we need to serialize A and B,  thus we need more information
     * about their length and an extra function. We can move these tasks to
     * to AbstractClawProblem. At the moment, we stick to the old method.
     * */

    static_assert(std::is_base_of<AbstractDomain<typename Domain_C::t>, Domain_C>::value,
		  "C not derived from AbstractDomain");
  }
  
  /* 
   * problem: to invoke the "is_equal(x, y)" method from the domain, we need an
   * object of type Domain_C, and we don't have any...
   */

  /* specification of the collision to find */


  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  
  /* assuming that f(x0) == f(x1) == y, is (x0, x1) an acceptable outcome? */
  bool is_good_pair(C_t const &z,  A_t const &x0,  A_t const &x1)
  { 
    return true;    // by default, yes.
  }
	
  /* embedding and randomization */

  /*
   * TODO: we need to formalize what we expect from this function.
   * It seems reasonable to ask that there must exist a subset D of C 
   * such that the embedding, when its domain is restricted to D, is injective. 
   */
  void send_C_to_A(const C_t &inp_C, A_t& out_A) const;
	
  /*
   * Family of permutations acting on the output domain.
   */
  void mix(const I &i, const C_t& x, C_t& y) const;   /* y <--- σ_i(x) */

  /* Generate a default permutation of C (e.g. the identity) */
  I mix_default() const; 

  /* Generate a new random permutation of C */
  I mix_sample(PRNG& rng) const; 
};
  
}

#endif

