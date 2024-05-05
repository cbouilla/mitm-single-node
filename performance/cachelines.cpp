/* This document to investigate  whether we should have three seperate arrays in
 * sender_buffers, or one large array of u8 that have triple of elements encoded.
 * compile:
 * g++ -O3 cachelines.cpp -o cachelines && time ./cachelines
 * According to godbolt compiler explorer, the compiler use vector instructions
 * for the three seperate arrays, while for the triple of elements single array
 * it doesn't use vector instructions!
 */


#include <iostream>
#include <algorithm>
#include <fstream>
#include <array>
#include <cstddef>
#include <cstdint>
#include<cstring>
#include <random>
#include <sys/types.h>
#include <vector>
#include <cstdlib>
#include <chrono>


double wtime() /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));


  return seconds;

}



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
    //gen64.seed(read_urandom<uint64_t>());
  };

  /* Initializes with a given a seed */
  PRNG(uint64_t seed){ gen64.seed(seed); }

  void update_seed() { gen64.seed(rd());  }

  /* Resets the internal state of */
  void set_seed(uint64_t s) { gen64.seed(s);}
  
  
  uint64_t rand(){
    return gen64();
    //return read_urandom<uint64_t>();
  };
  
};

int main()
{

  size_t N = 1LL<<25;


  int* a = new int[N];
  uint64_t* b = new uint64_t[N];
  uint16_t* c = new uint16_t[N];

  for (size_t i = 0; i < N; ++i){
    a[i] = i;
    b[i] = i*i+i;
    c[i] = (i^3)*i + 1;
  }
    
  PRNG prng{1};
  size_t idx;
  size_t nfetches = 1000;
  
  int aa = 0;
  uint64_t bb = 0;
  uint16_t cc = 0;

  double elapsed = 0;
  double t = 0;
  for (size_t k = 0; k < N; ++k){
    idx = prng.rand() % N;
    idx = std::min(idx, N - nfetches - 1); // not go out of bound
    
    t = wtime();
    for (size_t i = 0; i < nfetches; ++i) {
      aa += a[i] ^ i;
      bb += b[i] ^ i;
      cc += c[i] ^ i;
    }
    t = wtime() - t;
    elapsed += t;
    
  }



 

  std::cout << "three arrays:   elapsed = " << std::fixed << elapsed << "sec\n";
  std::cout << "ignore this line "<< aa << ", " << bb << ", " << cc << "\n";


  delete[] a;
  delete[] b;
  delete[] c;

  aa = 0;
  bb = 0;
  cc = 0;
  
  uint8_t* d = new uint8_t[N*(8 + 4 + 2)];
  prng.set_seed(1);



  for (size_t i = 0; i < N*14; ++i)
    d[i] = i;


  uint16_t* ptr_u16;
  uint64_t* ptr_u64;
  int* ptr_int;
  
  elapsed = 0;  
  for (size_t k = 0; k < N; ++k){
    idx = prng.rand() % N;
    idx = std::min(idx, N - (nfetches +1)*14); // not go out of bound
    
    t = wtime();
    /* Here we have an alignment issue to be investigated */ 
    for (size_t i = 0; i < nfetches; ++i) {
      ptr_u16 = (uint16_t*) &d[i*14];
      ptr_int = (int*) &d[i*14 + 2];
      ptr_u64 = (uint64_t*) &d[i*14 + 6];
      
      aa += *ptr_int ^ i;
      bb += *ptr_u64 ^ i;
      cc += *ptr_u16 ^ i;
    }
    t = wtime() - t;
    elapsed += t;
    
  }
  std::cout << "array(triples): elapsed = " << std::fixed << elapsed << "sec\n";
  std::cout << "ignore this line "<< aa << ", " << bb << ", " << cc << "\n";


  

}
