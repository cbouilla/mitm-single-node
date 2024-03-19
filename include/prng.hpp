#ifndef MITM_PRNG
#define MITM_PRNG

#include <fstream>
#include <array>
#include <cstddef>
#include <cstdint>
#include<cstring>
#include <random>

namespace mitm {

/*
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq
 */
class PRNG {
  /* source : https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine  */
  std::mersenne_twister_engine<std::uint_fast64_t, 64, 312, 156, 31,
			       0xb5026f5aa96619e9, 29,
			       0x5555555555555555, 17,
			       0x71d67fffeda60000, 37,
			       0xfff7eee000000000, 43,
			       6364136223846793005> gen64;
  std::random_device rd;

public:

  PRNG()
  { /* random seed */
    gen64.seed(rd());
  };

  void update_seed() { gen64.seed(rd());  }

  uint64_t rand(){ return gen64(); };
};


}

#endif
