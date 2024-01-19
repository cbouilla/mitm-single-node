#ifndef MITM_ABSTRACT_COLLISION
#define MITM_ABSTRACT_COLLISION
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
template<typename Domain_A,  typename Domain_C>
class AbstractCollisionProblem {
public:
  /* these lines have to be retyped again */
  using A_t = typename Domain_A::t;
  using C_t = typename Domain_C::t;

  static const int f_eq_g;
  
  AbstractCollisionProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
  }
  
  inline void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  inline void send_C_to_A(const C_t& inp_C, A_t& out_A) const;

  /* changes the behavior of the two above functions */
  void update_embedding(PRNG& rng); 
};

}

#endif
