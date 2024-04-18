#include "../mitm.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <vector>
#include "../include/AES.hpp"
#include "bits_lib.hpp"

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


#define CEIL(a, b) (((a) + (b)-1) / (b))

#define NBITS_A 20
#define NBITS_C 20

#define NBYTES_A CEIL(NBITS_A, 8)
#define NBYTES_C CEIL(NBITS_C, 8)
/* stop here. */



/*
 * repr is an encapsulation of whatever data used.
 */
struct SHA2_out_repr {
  // u8 data[NBYTES_C]; /* output */
  u8 data[64]; /* output */

  /* Constructor */
  SHA2_out_repr()  {
    for (size_t i = 0; i < NBYTES_C; ++i)
      data[i] = 0;

    /* These parts is always zero! */
    for (size_t i = NBYTES_A; i < 64; ++i)
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
    for (size_t i = 0; i < NBYTES_A; ++i)
      data[i] = other.data[i];
  }

  /* III. copy assignement  */
  SHA2_A_inp_repr& operator=(SHA2_A_inp_repr const& other) {
    for (size_t i = 0; i < NBYTES_A; ++i)
      data[i] = other.data[i];

    return *this;
  }


};






/* Implement << operator for SHA2_inp_repr and SHA2_out_repr */
std::ostream& operator<<(std::ostream& os, const SHA2_A_inp_repr& x)
{
  
  for (size_t i = 0; i < NBYTES_A;++i) {
    os << "0x" << std::setfill('0') << std::setw(2)
       <<  std::hex << static_cast<unsigned int>( x.data[i] ) << ", ";
  }
  return os;
}



std::ostream& operator<<(std::ostream& os, const SHA2_out_repr& x)
{
  for (size_t i = 0; i < NBYTES_C; ++i) {
    os << "0x" << std::setfill('0') << std::setw(2)
       <<  std::hex << static_cast<unsigned int>( x.data[i]) << ", ";
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
    constexpr u8 rem_bits = NBITS_C % 8;
    constexpr u8 mask = (1 << rem_bits) - 1;

    for(int i = 0; i < NBYTES_C; ++i )
      x.data[i] = prng.rand(); /* a bit overkill to call rand on a single byte! */

    /* remove execessive bits */
    if (rem_bits > 0){
      x.data[NBYTES_C - 1] = x.data[NBYTES_C - 1] & mask;
    }

  }

  
  inline
  bool is_equal(t const& x, t const& y) const
  {
    return ( bits_memcmp(x.data, y.data, NBITS_C) == 0);
  }

  inline
  void serialize(const t& in, u8* out) const
  {
    bits_memcpy(out, in.data, NBITS_C);
  }

  inline
  void unserialize(const u8* in, t& out) const
  {
    bits_memcpy(out.data, in, NBITS_C);
  }


  inline
  void copy(const t& in, t& out)
  {
    bits_memcpy(out.data, in.data, NBITS_C);
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
// AbstractClawProblem <I, A_t, C_t>
class SHA2_Problem
  : mitm::AbstractCollisionProblem<std::array<u8, 16>,
				   SHA2_A_inp_repr,
				   SHA2_OUT_DOMAIN>
{
public:
  SHA2_OUT_DOMAIN C; /* Output related functions */

  /* These types shorthand are essential ! */
  using I_t = std::array<u8, 16>; /* 1st tempalte type */
  using A_t = SHA2_A_inp_repr; /* 2nd template type */
  using C_t = SHA2_out_repr; /* 4th template type */
  using Dom_C = SHA2_OUT_DOMAIN;

  size_t nbits_A = NBITS_A;
  size_t nbits_C = NBITS_C;
  
  
  /****************************************************************************/
  // INIT
    SHA2_Problem()
  {
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, 255); /* a random byte value*/

    constexpr u8 rem_bits = NBITS_A % 8;
    constexpr u8 mask     =  (1 << rem_bits) - 1;
    
    /* Fill golden_input with uniform garbage, truly golden data. */
    for (int i = 0; i < NBYTES_A; ++i){
      golden_inp0.data[i] = distrib(gen);
      golden_inp1.data[i] = distrib(gen);
    }

    /* clear the rest of the bits */
    if (rem_bits > 0 ){
      golden_inp0.data[NBYTES_A - 1] = golden_inp0.data[NBYTES_A - 1] & mask;
      golden_inp1.data[NBYTES_A - 1] = golden_inp1.data[NBYTES_A - 1] & mask;
    }

    
      

    /* It's enough to compute golden_out for inp0, by def f(inp1) = golden_out */
    f(golden_inp0, golden_out);

    /* Check our constant work */
    C_t y0;
    C_t y1;
    f(golden_inp0, y0);
    f(golden_inp1, y1);
    
    std::cout << "\n========================================\n"
	      << "Does the golden pair collide? " << C.is_equal(y0, y1) << "\n"      
	      << "golden_inpA = " << golden_inp0 << "\n"
      	      << "golden_inpB = " << golden_inp1 << "\n"
	      << "golden_out  = " << golden_out  << "\n"
	      << "========================================\n";
  }
  /****************************************************************************/
  
  inline
  void f(const A_t& x, C_t& y) const
  {
    /* Artifically make golden_inp1 collides with golden_inp0 */
    if (is_equal_A(x, golden_inp1)){
      bits_memcpy(y.data, golden_out.data, NBITS_C);
      return;
    }

    /* Otherwise use normal sha2 */
    u32 state[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    sha256_process(state, x.data, 64);
    bits_memcpy(y.data, state, NBITS_C);
  }


  inline
  void send_C_to_A(C_t& inp_C, A_t& out_A) const
  {
    /* remove any junk in A */
    std::memset(out_A.data, 0, 64);
    bits_memcpy(out_A.data, inp_C.data, NBITS_C);
  }


  void mix(I_t const& key_static, C_t const& x, C_t& y) const
  {
    // std::cout << "before mix x = " << x << ", "
    // 	      << "y = " << y << "\n";

    char unsigned data[16] = {0};
    unsigned char key[16] = {0};

    std::memcpy(key, key_static.data(), 16);
    bits_memcpy(data, x.data, NBITS_C);

    Cipher::Aes<128> aes(key);
    aes.encrypt_block(data);
    bits_memcpy(y.data, data, NBITS_C);
  }

  
  I_t mix_default() const {return I_t {0}; }
  I_t mix_sample(mitm::PRNG& rng) const
  {
    I_t key_static{};
    /* get a "random" key */
    for (int k = 0; k <  16; ++k) 
      key_static[k] = rng.rand();

    return key_static;    
  }


  /* assuming that f(x0) == f(x1) == y, is (x0, x1) an acceptable outcome? */
  bool is_good_pair(C_t const &z,  A_t const &x0,  A_t const &x1)
  {
    /* TOO strong condition */
    // return (C.is_equal(golden_out, z)
    // 	    and is_equal_A(golden_inp0, x0)
    // 	    and is_equal_A(golden_inp1, x1));

    return C.is_equal(golden_out, z);
  }


public:
  /* A simple equality test for the type A_t, since it's not required by mitm */
  bool is_equal_A(A_t const& inp0, A_t const& inp1) const{
    return (bits_memcmp(inp0.data, inp1.data, NBITS_C) == 0); /* i.e. return not not_equal */
  }


public:
  C_t golden_out ; /* We look for point that equals this */
  A_t golden_inp0; /* We will edit f so that */
  A_t golden_inp1;
  C_t constant{};


};

/* This function enables us to compute how many bytes needed to have nslots
 * in the dictionary. The calculation method depends on how the dictionary is
 * implemented. It works for now on our simple dictionary.
 */
size_t calculate_nbytes_given(size_t nslots)
{
  /* sizeof(Value) + sizeof(Key) + sizeof(Chainlength)*/
  size_t triple_size = sizeof(SHA2_out_repr) + sizeof(u64) + sizeof(u64);

  return triple_size*nslots;
}


int main(int argc, char* argv[])
{
  
  SHA2_Problem Pb;

  if (argc == 2){
    /* ./sha2_collision_demo log2_nslots  */
    
    /*  Get the value of difficulty as a string */
    
    std::string log2_nslots_str = argv[1];
    size_t log2_nslots = std::stoull(log2_nslots_str);
    size_t nbytes_dict = calculate_nbytes_given(1LL<<log2_nslots);

    printf("nslots = 2^%lu\n", log2_nslots);

    /* Main work: */
    mitm::collision_search(Pb, nbytes_dict);

  } else if (argc == 3) {
    /* ./sha2_collision_demo log2_nslots difficulty */
    std::string log2_nslots_str = argv[1];
    std::string d_str = argv[2];
    int difficulty = std::stoi(d_str);
    size_t log2_nslots = std::stoull(log2_nslots_str);
    size_t nbytes_dict = calculate_nbytes_given(1LL<<log2_nslots);

    printf("nslots = 2^%lu, difficulty = %d\n", log2_nslots, difficulty);

    /* Main work */
    mitm::collision_search(Pb, nbytes_dict, difficulty);
  } else {
    mitm::collision_search(Pb);    
  }


}

