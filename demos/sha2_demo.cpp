#include "../mitm_sequential.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <array>

using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t ;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

#define NWORDS_DIGEST 8
#define WORD_SIZE 4 // bytes
#define NBYTES_DIGEST (NWORDS_DIGEST * WORD_SIZE)
/*
 * repr is an encapsulation of whatever data used.
 * It should have
 * I : Constructor
 * II: Destructor
 */
struct SHA2_out_repr {
  u32* state; /* output */

  /* Constructor */
  SHA2_out_repr() : state{new u32[NWORDS_DIGEST]} { }
  /* It could've been written as:
   * SHA2_out_repr(){ state = new u32[NBYTEST_DIGEST]; }
   */
  ~SHA2_out_repr()
  {
    delete[] state; 
  }
};

struct SHA2_inp_repr {
  u8* data; /* output */

  /* Constructor */
  SHA2_inp_repr() : data{new u8[64]} {}
  /* It could've been written as:
   * SHA2_inp_repr(){ state = new u8[64]; }
   */
  ~SHA2_inp_repr()
  {
    delete[] data;
  }
};



/*
 * Implement functions related to inp/out as specified in AbstractDomain
 */
class SHA2_OUT_DOMAIN : mitm::AbstractDomain<SHA2_out_repr>
{
  const static int length = NWORDS_DIGEST;
  const static size_t n_elements = (1LL<<32)<<length;
  /* todo: randomize */
  inline
  bool is_equal(SHA2_out_repr& x, SHA2_out_repr& y) const
  {
    return (0 == std::memcmp(x.state, y.state, NBYTES_DIGEST));
  }

  inline
  void serialize(const SHA2_out_repr& in, u8* out) const
  {
    std::memcpy(out, in.state, NBYTES_DIGEST);
  }

  inline
  void unserialize(const u8* in, SHA2_out_repr& out) const
  {
    std::memcpy(out.state, in, NBYTES_DIGEST);
  }


  inline
  void copy(const SHA2_out_repr& in, SHA2_out_repr& out)
  {
    std::memcpy(out.state, in.state, NBYTES_DIGEST);
  }

  inline
  int extract_1_bit(const SHA2_out_repr& inp) const
  {
    return (1&inp.state[0]);
  }

  inline
  u64 hash(const SHA2_out_repr& x) const
  {
    /* in case we are only extracting one word digest */
    constexpr size_t _2nd_idx = std::min(NWORDS_DIGEST, 1);

    return (static_cast<u64>(x.state[_2nd_idx]) << WORD_SIZE) | x.state[0];
  }

  
};

class SHA2_INP_DOMAIN : mitm::AbstractDomain<SHA2_inp_repr>
{
  const static int length = 64;
};







































