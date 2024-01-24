#ifndef MITM_ABSTRACT_COLLISION
#define MITM_ABSTRACT_COLLISION
#include "AbstractDomain.hpp"
#include <type_traits>
#include "include/prng.hpp"

namespace mitm {
/*
 * Provides the description of the types A, C and a function f : A -> C.
 *
 * The problem instance can contain extra data, e.g. if the goal consists
 * in finding H(prefix || x) == H(prefix || y) with x != y, then the Problem
 * could contain prefix.
 */
template<typename Domain_A,  typename Domain_C>
class AbstractCollisionProblem {
public:
  /* these lines have to be retyped again */
  using A_t = typename Domain_A::t;
  using C_t = typename Domain_C::t;
  
  AbstractCollisionProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
  }
  
  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  void send_C_to_A(const C_t& inp_C, A_t& out_A) const;

  /* assuming that f(x) == g(y) == z, is (x, y) an acceptable outcome? */
  bool good_collision(const A_t &x, const A_t &y, const C_t &z) const 
  { 
    /* 
     * problem: to invoke the "is_equal(x, y)" method from the domain, we need an
     * object of type Domain_A.
     */
    return (x != y); 
  }

  /* changes the behavior of the two above functions */
  /* CB: I am tempted to make the "index" of the function family explicit */
  void update_embedding(PRNG& rng); 
};

}

#endif
