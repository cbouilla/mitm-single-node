#include "../collision_engine.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <array>

/* We would like to call C function defined in `sha256.c` */
extern "C"{
 void sha256_process(uint32_t state[8], const uint8_t data[], uint32_t length);
 }




/* More aesthitically pleasing than uintXY_t */
using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t ;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

#define NWORDS_DIGEST 2
#define WORD_SIZE 4
#define NBYTES_A 64
#define NBYTES_B 64
#define NBYTES_DIGEST (NWORDS_DIGEST * WORD_SIZE)
/*
 * repr is an encapsulation of whatever data used.
 */
struct SHA2_out_repr {
  u32 state[NWORDS_DIGEST]; /* output */

  /* Constructor */
  SHA2_out_repr()  {
    for (size_t i = 0; i < NWORDS_DIGEST; ++i)
      state[i] = 0;
  }
};



struct SHA2_A_inp_repr {
  /* This is the active input */
  // u8 data[NBYTES_A]; /* output */
  /* but for convenience in writing f we make it as  following: */
  u8 data[64]; /* output */
  /* Constructor */
  SHA2_A_inp_repr()  {
    for (size_t i = 0; i < NBYTES_A; ++i)
      data[i] = 0;

    /* These parts is always zero! */
    for (size_t i = NBYTES_A; i < 64; ++i)
      data[i] = 0;

  }

};


struct SHA2_B_inp_repr {
  /* This is the active input */
  // u8 data[NBYTES_B]; /* output */
  /* but for convenience in writing f we make it as  following: */
  u8 data[64]; /* output */
  /* Constructor */
  SHA2_B_inp_repr()  {
    for (size_t i = 0; i < NBYTES_B; ++i)
      data[i] = 0;

    /* These parts is always zero! */
    for (size_t i = NBYTES_B; i < 64; ++i)
      data[i] = 0;
  }
};



/* Implement << operator for SHA2_inp_repr and SHA2_out_repr */
std::ostream& operator<<(std::ostream& os, const SHA2_A_inp_repr& x)
{
  
  for (size_t i = 0; i < NBYTES_A;++i) {
    os << "0x" << std::setfill('0') << std::setw(2) <<  std::hex << x.data[i] << ", ";
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, const SHA2_B_inp_repr& x)
{
  
  for (size_t i = 0; i < NBYTES_B;++i) {
    os << "0x" << std::setfill('0') << std::setw(2) <<  std::hex << x.data[i] << ", ";
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, const SHA2_out_repr& x)
{
  
  for (size_t i = 0; i < NWORDS_DIGEST; ++i) {
    os << "0x" << std::setfill('0') << std::setw(8) <<  std::hex << x.state[i] << ", ";
  }
  return os;
}




/*
 * Implement functions related to inp/out as specified in AbstractDomain
 */
class SHA2_OUT_DOMAIN : mitm::AbstractDomain<SHA2_out_repr>
{
public:
  // Not only a shortcut, it's necessry to avoid error when specializing a class
  // In template: 't' is a private member of 'mitm::AbstractDomain<SHA2_out_repr>'
  using t = SHA2_out_repr;
  
  const static int length = NBYTES_DIGEST;
  int a[NBYTES_DIGEST];
  
  const static size_t n_elements = (1LL<<length)*8;
  /* todo: randomize */
  inline
  void randomize(t& x, mitm::PRNG& prng) const
  {
    for(int i = 0; i<NWORDS_DIGEST; ++i )
      x.state[i] = prng.rand();
  }

  
  inline
  bool is_equal(t& x, t& y) const
  {
    return ( std::memcmp(x.state, y.state, NBYTES_DIGEST) == 0);
  }

  inline
  void serialize(const t& in, u8* out) const
  {
    std::memcpy(out, in.state, NBYTES_DIGEST);
  }

  inline
  void unserialize(const u8* in, t& out) const
  {
    std::memcpy(out.state, in, NBYTES_DIGEST);
  }


  inline
  void copy(const t& in, t& out)
  {
    std::memcpy(out.state, in.state, NBYTES_DIGEST);
  }

  inline
  int extract_1_bit(const t& inp) const
  {
    return (1&inp.state[0]);
  }

  inline
  u64 hash(const t& x) const
  {
    /* in case we are only extracting one word digest */
    constexpr size_t _2nd_idx = std::min(NWORDS_DIGEST - 1, 1);
    constexpr  size_t cond_shift = std::min(NWORDS_DIGEST - 1, 1);
    return (static_cast<u64>(x.state[_2nd_idx]) << 8*WORD_SIZE*cond_shift)
            | x.state[0];
  }

  
};


////////////////////////////////////////////////////////////////////////////////

class SHA2_Problem
  : mitm::AbstractClawProblem<uint64_t, SHA2_A_inp_repr, SHA2_B_inp_repr, SHA2_OUT_DOMAIN>
{
public:

  

  SHA2_OUT_DOMAIN C; /* Output related functions */

  /* These types shorthand are essential ! */
  using I_t = uint64_t; /* 1st tempalte type */
  using A_t = SHA2_A_inp_repr; /* 2nd template type */
  using B_t = SHA2_B_inp_repr; /* 3rd template type */
  using C_t = SHA2_out_repr; /* 4th template type */
  
  
  static const int f_eq_g = 0;
  
  inline
  void f(const A_t& x, C_t& y) const
  {
    u32 state[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    sha256_process(state, x.data, 64);
    std::memcpy(y.state, state, NBYTES_DIGEST);
  }


  inline
  void send_C_to_A(C_t& inp_C, A_t& out_A) const
  {
    /* remove any junk in A */
    std::memset(out_A.data, embedding_n, 64);
    std::memcpy(out_A.data, inp_C.state, NBYTES_DIGEST);
  }



  void mix(const I_t& i, const C_t& x, C_t& y) const {
    for (int j = 0; j<NWORDS_DIGEST; ++j)
      y.state[j] = x.state[j] ^ (i>>(j*32));
  }
  I_t mix_default() const {return 0;}
  I_t mix_sample(mitm::PRNG& rng) const {return rng.rand();}
  
private:
  /* Changes the extra bits in the input to embedding_n */
  int embedding_n = 0;
  
  
};


int main(int argc, char* argv[])
{
  SHA2_Problem Pb;
  mitm::collisoin_search(Pb);
}

