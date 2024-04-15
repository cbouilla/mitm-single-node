#include "../mitm.hpp"
#include "../include/AES.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <iomanip>
#include <array>
#include <string>

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
#define NBYTES_A 1 
#define NBYTES_B 1 
#define NBYTES_C 1 
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
    os << "0x" << std::setfill('0') << std::setw(2)
       <<  std::hex << static_cast<unsigned int>( x.data[i] ) << ", ";
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, const SHA2_B_inp_repr& x)
{
  
  for (size_t i = 0; i < NBYTES_B;++i) {
    os << "0x" << std::setfill('0') << std::setw(2)
       <<  std::hex << static_cast<unsigned int>(x.data[i]) << ", ";
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
    for(int i = 0; i < NBYTES_C; ++i )
      x.data[i] = prng.rand(); /* a bit overkill to call rand on a single byte! */
  }

  
  inline
  bool is_equal(t const& x, t const& y) const
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
// AbstractClawProblem <I, A_t, B_t, C_t>
class SHA2_Problem
  : mitm::AbstractClawProblem<std::array<u8, 16>,
			      SHA2_A_inp_repr,
			      SHA2_B_inp_repr,
			      SHA2_OUT_DOMAIN>
{
public:


  /* These types shorthand are essential ! */

  using A_t = SHA2_A_inp_repr; /* 2nd template type */
  using B_t = SHA2_B_inp_repr; /* 3rd template type */
  using C_t = SHA2_out_repr; /* 4th template type */
  using Dom_C = SHA2_OUT_DOMAIN;
  Dom_C C; /* Output related functions */  
  using I_t = std::array<u8, 16>;

  /****************************************************************************/
  // Initialization
  SHA2_Problem()
  {
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, 255); /* a random byte value*/

    
    /* Fill golden_input with uniform garbage, truly golden data. */
    for (int i = 0; i<NBYTES_A; ++i)
      golden_inpA.data[i] = distrib(gen);

    /* Fill golden_input with uniform garbage, truly golden data. */
    for (int i = 0; i<NBYTES_B; ++i){
      golden_inpB.data[i] = distrib(gen);
      constant.data[i] = 0;
    }


    /* get the constant, C,  that makes them collide */
    // f(golden_inpA) = g(golden_inpB) xor C
    f(golden_inpA, golden_out);
    g(golden_inpB, constant);

    
    /* our golden output is picked uniformly. Can `mitm` find a needle in
     * heystack? */
    for (int i = 0; i < NBYTES_C; ++i)
      constant.data[i] ^= golden_out.data[i];

    /* Check our constant work */
    C_t y0;
    C_t y1;
    f(golden_inpA, y0);
    g(golden_inpB, y1);
    
    std::cout << "\n========================================\n"
	      << "Does the golden pair collide? " << C.is_equal(y0, y1) << "\n"      
	      << "golden_inpA = " << golden_inpA << "\n"
      	      << "golden_inpB = " << golden_inpB << "\n"
	      << "golden_out  = " << golden_out  << "\n"
	      << "========================================\n";
  }
  /****************************************************************************/

  
  inline
  void f(const A_t& x, C_t& y) const
  {
    u32 data[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    sha256_process(data, x.data, 64);
    std::memcpy(y.data, data, NBYTES_C);

  }

  inline
  void g(const B_t& x, C_t& y) const
  {
    u32 data[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};

    u8* data_u8 = (u8*) data; /* a stupid trick to work on bytes */
    
    sha256_process(data, x.data, 64);
    /* Add our magical constant to get collision between:
     * golden_inpA and golden_inpB, i.e. f(golden_inpA) = g(golden_inpB)
     */
    for (int i = 0; i < NBYTES_C; ++i)
      data_u8[i] ^= constant.data[i];

    std::memcpy(y.data, data, NBYTES_C);
  }


  inline
  void send_C_to_A(C_t& inp_C, A_t& out_A) const
  {
    /* remove any junk in A */
    std::memset(out_A.data, 0, 64);
    /* Since domain and range may have different size: */
    std::memcpy(out_A.data, inp_C.data, std::min(NBYTES_A, NBYTES_C));
  }

  inline
  void send_C_to_B(C_t& inp_C, B_t& out_B) const
  {
    /* remove any junk in B */
    std::memset(out_B.data, 0, 64);
    /* Since domain and range may have different size: */
    std::memcpy(out_B.data, inp_C.data, std::min(NBYTES_B, NBYTES_C));
  }


  /* mix a value lives in C_t using i function scrambling. */
  void mix(I_t const& key_static, C_t const& x, C_t& y) const
  {
    // std::cout << "before mix x = " << x << ", "
    // 	      << "y = " << y << "\n";

    char unsigned data[16] = {0};
    unsigned char key[16] = {0};

    std::memcpy(key, key_static.data(), 16);
    std::memcpy(data, x.data, NBYTES_C);

    Cipher::Aes<128> aes(key);
    aes.encrypt_block(data);
    std::memcpy(y.data, data, NBYTES_C);
  
    // std::cout << "after mix x = " << x << ", "
    // 	      << "y = " << y << "\n";
  }
  
  I_t mix_default() const
  {
    return I_t {0};
  }
  
  I_t mix_sample(mitm::PRNG& rng) const
  {

    I_t key_static{};
    /* get a "random" key */
    for (int k = 0; k <  16; ++k) 
      key_static[k] = rng.rand();

    return key_static;    
  }

  bool is_good_pair(C_t const &z,  A_t const &x,  B_t const &y) const 
  {
    return (C.is_equal(golden_out, z)
	    and is_equal_A(golden_inpA, x)
	    and is_equal_B(golden_inpB, y));
  }

private:

public: /* they are public for debugging */
  C_t golden_out ; /* We look for point that equals this */
  A_t golden_inpA; /* We will edit f so that */
  B_t golden_inpB;
  C_t constant{};

  /* A simple equality test for the type A_t, since it's not required by mitm */
  bool is_equal_A(A_t const& inp0, A_t const& inp1) const{
    int not_equal = 0;
    for (int i = 0; i < NBYTES_A; ++i)
      not_equal += (inp0.data[i] != inp1.data[i]);

    return (not_equal == 0); /* i.e. return not not_equal */
  }

  /* A simple equality test for the type B_t, since it's not required by mitm */
  bool is_equal_B(B_t const& inp0, B_t const& inp1) const{
    int not_equal = 0;
    for (int i = 0; i < NBYTES_B; ++i)
      not_equal += (inp0.data[i] != inp1.data[i]);

    return (not_equal == 0); /* i.e. return not not_equal */
  }

  /* print if the value in C equals the golden input in A or B */
  void test_C_eq_golden_inp(C_t&)
  {
    
  }

};


int main(int argc, char* argv[])
{

  mitm::nbytes_A = NBYTES_A;
  mitm::nbytes_B = NBYTES_B;
  SHA2_Problem Pb;
  
  std::cout << "sha2-claw demo! |inp_A| = " << NBYTES_A << "bytes, "
	    << "|inp_B| = " << NBYTES_B << "bytes, "
	    << "|out| = " << NBYTES_C << "bytes\n";


  if (argc < 2) /* we have not provided difficulty */
    mitm::claw_search(Pb);
  if (argc == 2){
    /*  Get the value of difficulty as a string */
    std::string d_str = argv[1];
    int difficulty = std::stoi(d_str);
    mitm::claw_search(Pb, difficulty);
  }
    
}

