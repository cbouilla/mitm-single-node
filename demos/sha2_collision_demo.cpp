#include "../mitm.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>

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

/* Edit the next *three* line to define a new demo! */
#define NBYTES_A 2 /* A's input length <= 64 bytes  */
#define NBYTES_C 2 /* output length <= 32 bytes */
/* stop here. */



/*
 * repr is an encapsulation of whatever data used.
 */
struct SHA2_out_repr {
  u8 data[NBYTES_C]; /* output */

  /* Constructor */
  SHA2_out_repr()  {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = 0;
  }

  /* I. destructor */
  ~SHA2_out_repr(){}; /* no dynamic memory to be removed */
  
  /* II. copy constructor */
  SHA2_out_repr(SHA2_out_repr const& other) {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = other.data[i];
  }

  /* III. copy assignement  */
  SHA2_out_repr& operator=(SHA2_out_repr const& other) {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = other.data[i];

    return *this;
  }


  /****************************************************************************/
  // For the naive algorithm
  bool operator<(SHA2_out_repr const& other) const {
    return (std::memcmp(data, other.data, NBYTES_C) < 0);
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

    /* II. copy constructor */
  SHA2_A_inp_repr(SHA2_A_inp_repr const& other) {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = other.data[i];
  }

  /* III. copy assignement  */
  SHA2_A_inp_repr& operator=(SHA2_A_inp_repr const& other) {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = other.data[i];

    return *this;
  }

  /****************************************************************************/
  // For the naive algorithm
  bool operator<(SHA2_A_inp_repr const& other) const {
    return (std::memcmp(data, other.data, NBYTES_A) < 0);
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



std::ostream& operator<<(std::ostream& os, const SHA2_out_repr& x)
{
  os << "0x" ;
  for (size_t i = 0; i < NBYTES_C; ++i) {
    os << std::setfill('0') << std::setw(1) <<  std::hex << x.data[i] << ", ";
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
  
  const static int length = NBYTES_C;
  int a[NBYTES_C];
  
  const static size_t n_elements = (1LL<<length)*8;
  /* todo: randomize */
  inline
  void randomize(t& x, mitm::PRNG& prng) const
  {
    for(int i = 0; i < NBYTES_C; ++i )
      x.data[i] = prng.rand(); /* a bit overkill to call rand on a single byte! */
  }

  
  inline
  bool is_equal(t& x, t& y) const
  {
    return ( std::memcmp(x.data, y.data, NBYTES_C) == 0);
  }

  inline
  void serialize(const t& in, u8* out) const
  {
    std::memcpy(out, in.data, NBYTES_C);
  }

  inline
  void unserialize(const u8* in, t& out) const
  {
    std::memcpy(out.data, in, NBYTES_C);
  }


  inline
  void copy(const t& in, t& out)
  {
    std::memcpy(out.data, in.data, NBYTES_C);
  }

  inline
  int extract_1_bit(const t& inp) const
  {
    return (1&inp.data[0]);
  }

  inline
  u64 hash(const t& x) const
  {
    /* Take into account all bytes of the output */
    /* each bit k in position 64*r + l, contributes to the l poisition in the digest */
    u64 digest = 0;
    for (int i = 0; i < NBYTES_C; ++i)
      digest ^= ((u64) x.data[i])<<((8*i)%64);

    return digest;
  }


  
};


////////////////////////////////////////////////////////////////////////////////

class SHA2_Problem
  : mitm::AbstractCollisionProblem<uint64_t, SHA2_A_inp_repr, SHA2_OUT_DOMAIN>
{
public:
  SHA2_OUT_DOMAIN C; /* Output related functions */

  /* These types shorthand are essential ! */
  using I_t = uint64_t; /* 1st tempalte type */
  using A_t = SHA2_A_inp_repr; /* 2nd template type */
  using C_t = SHA2_out_repr; /* 4th template type */
  using Dom_C = SHA2_OUT_DOMAIN;
  
  static const int f_eq_g = 0;
  
  static inline
  void f(const A_t& x, C_t& y)
  {
    u32 state[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    sha256_process(state, x.data, 64);
    std::memcpy(y.data, state, NBYTES_C);
  }


  inline
  void send_C_to_A(C_t& inp_C, A_t& out_A) const
  {
    /* remove any junk in A */
    std::memset(out_A.data, embedding_n, 64);
    std::memcpy(out_A.data, inp_C.data, NBYTES_C);
  }



  void mix(const I_t& i, const C_t& x, C_t& y) const {
    for (int j = 0; j < NBYTES_C; ++j)
      y.data[j] = x.data[j] ^ (i>>(j*8));
  }
  I_t mix_default() const {return 0;}
  I_t mix_sample(mitm::PRNG& rng) const {return rng.rand();}


  void collect_all_collisions_naive()
  { /* Get all collisions by the naive method. */
    all_collisions = mitm::naive_collisoin_search(*this,
						  [](SHA2_A_inp_repr& x, u64 i)
						  {*((u64*) &x.data[0]) = i; },
						  (1LL<<(NBYTES_A*8))// |A|
						  );
    all_collisions_collected = true;
  }

private:
  /* Changes the extra bits in the input to embedding_n */
  int embedding_n = 0;
  /*----------------------------------------------------------------------------
   * In this demo, we only say the triple (z, x0, x1) is a golden collision
   * after we observe all the possible triples that makes collision i.e.
   * z = f(x0) = f(x1). We know we have exhausted all triples, since we computed
   * them using the naïve method. You can define is_good_collision as you wish!
   *
   ******************************************************************************/
  bool all_collisions_collected = false;
  std::vector<std::pair<SHA2_out_repr, SHA2_A_inp_repr>> all_collisions{};
  

};


int main()
{
  
  SHA2_Problem Pb;
  Pb.collect_all_collisions_naive();
  mitm::collision_search(Pb);
}

