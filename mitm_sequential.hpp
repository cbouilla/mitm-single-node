#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <array>

#include <type_traits>
#include<cstring>
#include <ios>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <assert.h>
#include <tuple>
#include <utility>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <string>

#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "AbstractCollisionProblem.hpp"
#include "claw_engine.hpp"
#include "collision_engine.hpp"
#include "engine.hpp"
#include "dict.hpp"
#include "include/memory.hpp"
#include "include/timing.hpp"
#include "include/prng.hpp"

namespace mitm {






template<typename Problem>
std::pair<typename Problem::C_t, typename Problem::C_t> collision(Problem& Pb) 
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;
  PRNG rng_urandom;
  
  /* Sanity Test: */
  std::cout << "unserial(serial(.)) =?= id(.) : "
	    << is_serialize_inverse_of_unserialize<Problem>(Pb, rng_urandom)
	    << "\n";


  // --------------------------------- INIT -----------------------------------/
  size_t n_bytes = 0.5*get_available_memory();
  
  std::cout << "Going to use "
	    << std::dec << n_bytes << " bytes = 2^"<< std::log2(n_bytes)
	    << " bytes for dictionary!\n";
  
  Dict<u64, C_t, Problem> dict{n_bytes}; /* create a dictionary */
  std::cout << "Initialized a dict with " << dict.n_slots << " slots = 2^"
	    << std::log2(dict.n_slots) << " slots\n";


  // -----------------------------------------------------------------------------/
  // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
  int difficulty = 9; // difficulty;
  /* inp/out variables are used as input and output to save one 1 copy */


  /***************************************************************/
  /* when generating a distinguished point we have:              */
  /*  1)   inp0           =f/g=> out0                            */
  /*  2)  (inp1 := out0)  =f/g=> out1                            */
  /*  3)  (inp2 := out1)  =f/g=> out2                            */
  /*            ...                                              */
  /* m+1) (inp_m := out_m) =f/g=> out_m                          */
  /* A distinguished point found at step `m+1`                   */
  /* "We would like to preserve inp0 at the end of calculation." */
  /* In order to save ourselves from copying in each step.       */
  /***************************************************************/

  C_t pre_inp0{}; /* We need to save this value, see above. */

  /* 1st set of buffers: Related to input0 as a starting point */
  /* either tmp0 or  output0 */
  C_t inp0_or_out0_buffer0{};
  C_t inp0_or_out0_buffer1{};

  C_t* inp0_pt    = &inp0_or_out0_buffer0;
  /* Always points to the region that contains the output */
  C_t* out0_pt = &inp0_or_out0_buffer1;
  u64  out0_digest = 0; /* hashed value of the output0 */
  /* Recall: Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  u8 out0_bytes[Pb.C.length];

  /* 2nd set of buffers: Related to input1 as a starting point */
  /* When we potentially find a collision, we need 2 buffers for (inp1, out1) */
  /* We will use the initial value of inp1 once, thus we don't need to gaurd  */
  /* in an another variable untouched */
  C_t inp1_or_out1_buffer0{};
  C_t inp1_or_out1_buffer1{};
  C_t* inp1_pt    = &inp1_or_out1_buffer0;
  /* Always points to the region that contains the output */
  C_t* out1_pt = &inp1_or_out1_buffer1;

  /* Use these variables to print the full collision */
  A_t inp_A{};
  B_t inp_B{};
  u8 inp_A_serial[Pb.A.length];
  u8 inp_B_serial[Pb.B.length];

  /* Store the results of collisions here */
  /* a:A_t -f-> x <-g- b:B_t */ 
  std::vector< std::pair<A_t, B_t> >  collisions_container{};



  /**************************** Collisions counters ***************************/
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  size_t chain_length1 = 0;
  
  bool is_collision_found = false;
  size_t n_collisions = 0;
  size_t n_needed_collisions = 1LL<<20;
  
  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  size_t n_robinhoods = 0;
  

  /*********************************************************
   * surprise we're going actually use
   * void iterate_once(typename Problem::C_t* inp_pt,
   *		  typename Problem::C_t* out_pt,
   *		  Iterate_F<Problem>& F,
   *		  Iterate_G<Problem>& G)
   *
   *********************************************************/

  bool found_dist = false;
  size_t n_distinguished_points = 0;
  constexpr size_t interval = (1LL<<15);
  double collision_timer = wtime();
  
  std::cout << "about to enter a while loop\n";
  /*------------------- Generate Distinguished Points ------------------------*/
  // while (n_collisions < 1){
  while (n_collisions < n_needed_collisions){

    /* These simulations show that if 10w distinguished points are generated
     * for each version of the function, and theta = 2.25sqrt(w/n) then ...
     */
    /* update F and G by changing `send_C_to_A` and `send_C_to_B` */    
    Pb.update_embedding(rng_urandom);
    dict.flush();

    for (size_t n_dist_points = 0;
	 n_dist_points < 10*(dict.n_slots);
	 ++n_dist_points)
      {
      is_collision_found = false;
      /* fill the input with a fresh random value. */
      Pb.C.randomize(pre_inp0, rng_urandom);
      chain_length0 = 0; /*  */

      found_dist = generate_dist_point<Problem>(pre_inp0,
						inp0_pt,
						out0_pt,
						chain_length0,
						difficulty,
						Pb);
      out0_digest = Pb.C.hash(*out0_pt);
      ++n_distinguished_points;

      print_interval_time(n_distinguished_points, interval);
      
      if (not found_dist) [[unlikely]]
	continue; /* skip all calculation below and try again  */
      
      ++n_dist_points;

      is_collision_found = dict.pop_insert(out0_digest, /* key */
					   pre_inp0, /* value  */
					   chain_length0,
					   *inp1_pt, /* save popped element here */
					   chain_length1, /* of popped input */
					   Pb);
      
      if (is_collision_found) [[unlikely]]{
	++n_collisions;

	/* Move this code to print collision information */
        std::cout << "\nA collision is found\n"
		  << "It took " << (wtime() - collision_timer) << " sec\n"
		  << "inp0 (starting point) = " << pre_inp0 << "\n"
		  << "digest0 = 0x" << out0_digest << "\n"
		  << "chain length0 = " << chain_length0 << "\n"
		  << "inp1 (starting point) = " << *inp1_pt << "\n"
	  	  << "chain length1 = " << std::dec << chain_length1 << "\n"
		  << "-------\n";

	collision_timer = wtime();
	/* respect the rule that inp0 doesn't have pointers dancing around it */
	Pb.C.copy(pre_inp0, *inp0_pt); /* (*tmp0_ptO) holds the input value  */

	/* i.e. when f =/= g then one of the inputs has to  correspond to A
	 * and the other has to correspond to B. The order doesn't matter.
	 * todo: it should also neglect robinhood.
	 */
	bool is_potential_coll = treat_collision<Problem>(inp0_pt,
							  out0_pt,
							  chain_length0,
							  inp1_pt,
							  out1_pt,
							  chain_length1,
							  collisions_container,
							  Pb);

	/* move print_collision_info here */

	bool real_collision = Pb.C.is_equal(*out0_pt, *out1_pt);
	bool is_robinhood = Pb.C.is_equal(*inp0_pt, *inp1_pt);
	n_robinhoods += is_robinhood;

	if (is_potential_coll){
	  /* remove this printing */
	  std::cout << "After treating collision\n"
		    << "inp0 = " << *inp0_pt << "\n"
		    << "out0 = " << *out0_pt << "\n"
		    << "inp1 = " << *inp1_pt << "\n"
		    << "out1 = " << *out1_pt << "\n"
		    << "out0 == out1? " << real_collision  << "\n"
		    << "diges0 == digest1? "
		    << "robinhood? " << is_robinhood  << "\n"
		    << "#collisions = " << std::dec << n_collisions << "\n"
		    << "#robinhood = "  << std::dec << n_robinhoods << "\n"
		    << "\n";


	  /* Get the complete inputs as they live in A and B */
	  std::cout << "container length " << collisions_container.size() << "\n"
		    << "is a good collisision? " << is_potential_coll << "\n";
	  Pb.A.serialize(collisions_container.back().first, inp_A_serial);
	  Pb.B.serialize(collisions_container.back().second, inp_B_serial);

	  printf("inp_A = {");
	  for(size_t j = 0; j < Pb.A.length; ++j)
	    printf("0x%02x, ", inp_A_serial[j]);
	  puts("};");

	  printf("inp_B = {");
	  for(size_t j = 0; j < Pb.B.length; ++j)
	    printf("%02x, ", inp_B_serial[j]);
	  puts("};\n________________________________________\n");
	}
      }
    }
  }
  /* end of work */
  return std::pair<C_t, C_t>(*inp0_pt, *inp1_pt); // todo wrong values
}


template <typename Problem>
void find_claw(Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;

  C_t pre_inp0{}; /* We need to save this value, see above. */

  /* 1st set of buffers: Related to input0 as a starting point */
  /* either tmp0 or  output0 */
  C_t inp0_or_out0_buffer0{};
  C_t inp0_or_out0_buffer1{};

  C_t* inp0_pt    = &inp0_or_out0_buffer0;
  /* Always points to the region that contains the output */
  C_t* out0_pt = &inp0_or_out0_buffer1;
  u64  out0_digest = 0; /* hashed value of the output0 */
  /* Recall: Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  u8 out0_bytes[Pb.C.length];

  /* 2nd set of buffers: Related to input1 as a starting point */
  /* When we potentially find a collision, we need 2 buffers for (inp1, out1) */
  /* We will use the initial value of inp1 once, thus we don't need to gaurd  */
  /* in an another variable untouched */
  C_t inp1_or_out1_buffer0{};
  C_t inp1_or_out1_buffer1{};
  C_t* inp1_pt    = &inp1_or_out1_buffer0;
  /* Always points to the region that contains the output */
  C_t* out1_pt = &inp1_or_out1_buffer1;

  /* Use these variables to print the full collision */
  A_t inp_A{};
  B_t inp_B{};
  u8 inp_A_serial[Pb.A.length];
  u8 inp_B_serial[Pb.B.length];

  /* Store the results of collisions here */
  /* a:A_t -f-> x <-g- b:B_t */ 
  std::vector< std::pair<A_t, B_t> >  collisions_container{};

}

// to use parallel sort sort
// install tbb lib
// sudo apt install libtbb-dev
// compiling
// g++ -ggdb -O0 -std=c++17 demos/speck32_demo.cpp -o speck32_demo -ltbb
// g++ -flto -O3 -std=c++17  -fopenmp demos/speck32_demo.cpp -o speck32_demo -ltbb

  
}

#endif
