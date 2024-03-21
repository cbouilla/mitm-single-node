#ifndef MITM_ABSTRACT_DOMAIN
#define MITM_ABSTRACT_DOMAIN
#include <cstdint>
#include <cstddef>

namespace mitm
{
using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t ;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

/******************************************************************************/
/* Document for standard implementation                                       */
/******************************************************************************/
/*
 * A "domain" extends some type to provide the extra functions we need.
 * An instance of the domain can contain extra information
 * E.g. for small integers mod p, "repr" could be u64 while the domain actually contains p
 * E.g. for points on an elliptic curve, "repr" could be a pair of integers mod p, while the
 *      domain would contain the equation of the curve, etc.
 */
template<class repr>           /* repr must support comparisons, and assignment */
class AbstractDomain {
public:
  int length;  /* nbytes needed to encode an element */
  size_t n_elements; /* how many elements in the domain */
  using t = repr;            /* t is the machine representation of elements of the domain */

  template<typename PRNG>
  void randomize(t &x, PRNG &p) const;           /* set x to a random value */

  bool is_equal(const t &x, const t &y) const;

  void serialize(const t &x, u8 *out) const;   /* write this to out */
  void unserialize(const u8 *in, t &x) const;        /* read this from in */
  void copy(const t& inp, t& out) const; /* deepcopy inp to out */


  u64 hash(const t &x) const;                /* return some bits from this */
  u64 hash_extra(const t &x) const ;         /* return more bits from this */
};


}

#endif
