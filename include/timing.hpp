#ifndef MITM_ACCESSSOIRES
#define MITM_ACCESSSOIRES
#include <chrono>
#include <fstream>
#include <array>
#include <cstddef>
#include <cstdint>
#include<cstring>
#include <ios>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <sys/types.h>


namespace mitm
{

  

static /* since it's defined in a header file */
double wtime() /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));


  return seconds;

}

 /*
  * Print how long does it took to do n iterations.
  */
inline static void print_interval_time(size_t n)
{
  static double previous_time = -1;
  double current_time = wtime();
  double elapsed = current_time - previous_time;
  
  if (previous_time == -1){
    previous_time = wtime();
    return; /* 1st measurement */
  }


  printf("\r%lu=2^%0.2f iter took:"
	 " %0.2f sec, i.e. %0.2f = 2^%0.2f iter/sec",
	 n, std::log2(n),
	 elapsed, n/elapsed, std::log2(n/elapsed) );
  fflush(stdout);
  previous_time = wtime();
  
  
}


}    
#endif
