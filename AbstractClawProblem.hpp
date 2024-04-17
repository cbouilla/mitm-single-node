#ifndef MITM_ABSTRACT_CLAW
#define MITM_ABSTRACT_CLAW

#include "AbstractDomain.hpp"
#include <type_traits>
#include "include/prng.hpp"

namespace mitm {
  /*
   * Provides the description of the types A, B and C, as well as 
   * two functions f : A -> C and g : B -> C and an optionnal
   * predicate P.
   *
   * The problem instance can contain extra data.
   * E.g. In the attack on double-encryption, where the goal is to find
   *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
   */
template<typename I, typename A, typename B, typename Domain_C>
class AbstractClawProblem {
public:
  /* these lines have to be retyped again */
  using I_t = I;
  using A_t = A;
  using B_t = B;
  using C_t = typename Domain_C::t;
  using Dom_C = Domain_C;
  Dom_C C;
  /* Used in collecting stats */
  size_t nbits_A;
  size_t nbits_B;
  size_t nbits_C;

  AbstractClawProblem() {
    static_assert(std::is_base_of<AbstractDomain<typename Domain_C::t>, Domain_C>::value,
		  "C not derived from AbstractDomain");
  }
  
  /* specification of the collision to find */

  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  void g(const B_t &x, C_t &y) const;  /* y <--- g(x) */

  /* assuming that f(x) == g(y) == z, is (x, y) an acceptable outcome? */
  bool is_good_pair(C_t const &z,  A_t const &x,  B_t const &y)
  { 
    return true;    // by default, yes.
  }

  /* embedding and randomization */

  /*
   * TODO: we need to formalize what we expect from these two functions.
   * It seems reasonable to ask that there must exist a subset D of C 
   * such that these two, when their domain is restricted to D, are injective. 
   */
  void send_C_to_A(const C_t &inp_C, A_t &out_A) const;
  void send_C_to_B(const C_t &inp_C, B_t &out_B) const;

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
