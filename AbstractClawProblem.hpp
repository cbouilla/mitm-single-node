#ifndef MITM_ABSTRACT_CLAW
#define MITM_ABSTRACT_CLAW

#include "AbstractDomain.hpp"
#include <type_traits>
#include "include/prng.hpp"

namespace mitm {
/*
 * Provides the description of the type A and a function f : A -> A.
 *
 * The problem instance can contain extra data.
 * E.g. In the attack on double-encryption, where the goal is to find
 *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
 */
template<typename Domain_A, typename Domain_B, typename Domain_C>
class AbstractClawProblem {
public:
  /* these lines have to be retyped again */
  using C_t = typename Domain_C::t;
  using A_t = typename Domain_A::t;
  using B_t = typename Domain_B::t;
  
  AbstractClawProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
    static_assert(std::is_base_of<AbstractDomain<typename Domain_B::t>, Domain_B>::value,
      "B not derived from AbstractDomain");
  }
  
  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  void g(const B_t &x, C_t &y) const;  /* y <--- g(x) */
  void send_C_to_A(const C_t& inp_C, A_t& out_A) const;
  void send_C_to_B(const C_t& inp_C, B_t& out_B) const;

  /* assuming that f(x) == g(y) == z, is (x, y) an acceptable outcome? */
  bool good_collision(const A_t& x, const B_t& y, const z_t& z) const 
  { 
    return true; 
  }

  /* changes the behavior of the two above functions */
  void update_embedding(PRNG& rng); 
};



}

#endif
