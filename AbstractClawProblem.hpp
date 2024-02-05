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
template<typename I, typename Domain_A, typename Domain_B, typename Domain_C>
class AbstractClawProblem {
public:
  /* these lines have to be retyped again */
  using A_t = typename Domain_A::t;
  using B_t = typename Domain_B::t;
  using C_t = typename Domain_C::t;
  using I_t = I;
  
  AbstractClawProblem() {
    /* At the end, we need to serialize A and B,  thus we need more information
     * about their length and an extra function. We can move these tasks to
     * to AbstractClawProblem. At the moment, we stick to the old method.
     * */
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
    static_assert(std::is_base_of<AbstractDomain<typename Domain_B::t>, Domain_B>::value,
		  "B not derived from AbstractDomain");
    static_assert(std::is_base_of<AbstractDomain<typename Domain_C::t>, Domain_C>::value,
		  "C not derived from AbstractDomain");

  }
  
  /* specification of the collision to find */

  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  void g(const B_t &x, C_t &y) const;  /* y <--- g(x) */
        
  /* assuming that f(x) == g(y) == z, is (x, y) an acceptable outcome? */
  bool good_collision(const A_t &x, const B_t &y, const C_t &z) const 
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
