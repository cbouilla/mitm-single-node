/// Primitive adaptation of speck blockcipher using 16-bit words.
/// speck32* : {0, 1}^32 x {0, 1}^32 -> {0, 1}^32, the official documentation
/// specifies different length for k, we are using shorter key for simplicity
/// since security is not a concern.

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdint.h>
#include <iostream>
#include <array>
#include "../mitm_sequential.hpp"
#include "speck32.hpp"
#include <iomanip>

using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t ;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;



using Speck32_t = typename std::array<uint16_t , 2>;


std::ostream& operator<<(std::ostream& os, const Speck32_t& inp )
{
  os << "[0x"
     << std::hex << std::setfill('0') << std::setw(4) << inp[0]
     << ", 0x"
     << std::hex << std::setfill('0') << std::setw(4) << inp[1] << "]";
  return os;
}

class SPECK_DOMAIN : mitm::AbstractDomain< Speck32_t >{
  
public:
  using t = Speck32_t;
  static const int length = 4;
  static const size_t n_elements = (1LL<<32);
  
  inline bool is_equal(const t& x, const t& y) const
  { return x == y;  }

  void randomize(t& x, mitm::PRNG& prng) const
  {
    x[0] ^= prng.rand();
    x[1] ^= (prng.rand())>>16;
  }

  

  void serialize(const t& x, uint8_t * out) const
  {
    out[0] = x[0];
    out[1] = x[0]>>8;
    out[2] = x[1];
    out[3] = x[1]>>8;
  }

  void serialize(const t& x, std::array<uint8_t, length>& out) const
  {
    out[0] = x[0];
    out[1] = x[0]>>8;
    out[2] = x[1];
    out[3] = x[1]>>8;
  }

  void unserialize(const uint8_t* in, t& out) const
  {
    out[0] = in[0] | ((uint16_t ) in[1])<<8;
    out[1] = in[2] | ((uint16_t ) in[3])<<8;
  }

  void unserialize(const std::array<uint8_t, length>& in, t& out) const
  {
    out[0] = in[0] | ((uint16_t ) in[1])<<8;
    out[1] = in[2] | ((uint16_t ) in[3])<<8;
  }

  
  inline auto extract_1_bit(const t& inp) const
  { return inp[0]&1;  }


  inline uint32_t hash(const t& inp) const
  { return inp[0] | ((uint32_t) inp[1])<<16; }
  

  inline void copy(const t& inp, t& out) const
  {
      out[0] = inp[0];
      out[1] = inp[1];
  }

  inline void next(t& inp) const
  {
    uint32_t n = inp[0] | ((uint32_t) inp[1])<<16;
    ++n;
    inp[0] = n;
    inp[1] = n>>16;
  }

  inline void ith_elm(t& inp_out, size_t i) const
  {
    inp_out[0] = i;
    inp_out[1] = i>>16;
  }

};

class Problem : mitm::AbstractProblem<SPECK_DOMAIN, SPECK_DOMAIN, SPECK_DOMAIN >{
private:
  u64 version_send_C_to_A = 0;
  u64 version_send_C_to_B = 0;

public:
  using Dom_A = SPECK_DOMAIN;
  using Dom_B = SPECK_DOMAIN;
  using Dom_C = SPECK_DOMAIN;
  SPECK_DOMAIN A;
  SPECK_DOMAIN B;
  SPECK_DOMAIN C;
  using A_t = typename SPECK_DOMAIN::t;
  using B_t = typename SPECK_DOMAIN::t;
  using C_t = typename SPECK_DOMAIN::t;

  void f(const A_t &x, C_t &y) const
  {
    static std::array<uint16_t, 2> inpt_text{0, 0};
    encrypt(inpt_text, y, x);
  }

  void g(const B_t &x, C_t &y)  const
  {
    static std::array<uint16_t, 2> inp_ciphertext{0, 0};
    decrypt(inp_ciphertext, y, x);
  }

  void send_C_to_A(C_t& inp_C, A_t& out_A) const
  {
    out_A = inp_C;
    out_A[0] ^= version_send_C_to_A;
    out_A[0] ^= (version_send_C_to_A>>16);
  }
  void send_C_to_B(C_t& inp_C, B_t& out_B) const
  {
    out_B = inp_C;
    out_B[0] ^= version_send_C_to_B;
    out_B[0] ^= (version_send_C_to_B>>16);
  }

  void update_embedding(mitm::PRNG& rng) {
    /*
     * Change `send_C_to_A` and `send_C_to_B` functions
     */
    version_send_C_to_A = rng.rand();
    version_send_C_to_B = rng.rand();
  }
  
};



int main(int argc, char* argv[])
{
  Problem Pb;
  mitm::collision<Problem>(Pb);
}

