#ifndef MITM_MPI_RECEIVER
#define MITM_MPI_RECEIVER

#include <cstddef>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <mpi.h>

#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "AbstractCollisionProblem.hpp"
#include "include/prng.hpp"
#include "include/memory.hpp"
#include "include/counters.hpp"
#include "dict.hpp"
#include "engine.hpp"
#include "mpi_common.hpp"


namespace mitm {
/* One round of computation of receiver  */
template <typename Problem, typename... Types> bool receiver_round() {}



template<typename Problem,  typename... Types>
void reciever(Problem& Pb,
	      MITM_MPI_data& my_info,
	      int difficulty,
	      typename Problem::C_t& inp_St, // Startign point in chain
	      typename Problem::C_t* inp0_pt,
	      typename Problem::C_t* inp1_pt,
	      typename Problem::C_t* out0_pt,
	      typename Problem::C_t* out1_pt,
	      typename Problem::C_t& inp_mixed,
	      typename Problem::A_t& inp0A,
	      typename Problem::A_t& inp1A,
	      Types... args) /* two extra arguments if claw problem */
{
  using C_t = typename Problem::C_t;
  size_t const buf_size = my_info.nelements_buffer * Pb.C.size;
  // ======================================================================== //
  //                          STEP 0:  INTITALIZATION                         //
  // ------------------------------------------------------------------------ //

  
  /* ------------------------------- DICT INIT ------------------------------ */
  size_t nbytes_memory = get_nbytes_per_receivers(my_info);
  
  Dict<u64, C_t, Problem> dict{nbytes_memory};
  printf("Recv %d initialized a dict with %llu slots = 2^%0.2f slots\n",
	 my_info.my_rank_intercomm, dict.n_slots, std::log2(dict.n_slots));

  /* ----------------------------- receive buffer ---------------------------- */
  u8* rcv_buf  = new u8[buf_size];
  u8* work_buf = new u8[buf_size];

  
  /* ---------------- data inserted/extracted from dictionarry -------------- */
  // u64 out0_digest = 0;
  // size_t chain_length1 = 0;
  // // C_t* out1_pt,
  
  // /*--------------------------- Collisions counters --------------------------*/
  // /* How many steps does it take to get a distinguished point from  an input */
  // size_t chain_length0 = 0;
  // Counters ctr{}; /* various non-essential counters */

  // /* ----------------------------- Query Results ---------------------------- */
  // bool found_a_collision = false;
  
  // /* We should have ration 1/3 real collisions and 2/3 false collisions */
  // bool found_dist = false;

  // /* we found a pair of inputs that lead to the golden collisoin or golden claw! */
  // bool found_golden_pair = false;



  // ======================================================================== //
  //                                 STEP 1                                   //
  //      Initialization, and  agreeing on mix_function and walk function     //
  // ------------------------------------------------------------------------ //
  u64 mixer_seed;
  u64 byte_hasher_seed;

  seed_agreement(my_info, mixer_seed, byte_hasher_seed);

  /* Now, mixing function and byte_hasher are the same among all processes */
  PRNG prng_mix{mixer_seed}; /* for mixing function */
  PearsonHash byte_hasher{byte_hasher_seed}; /* todo update PearsonHash to accept a seed */
  /* variable to generate families of functions f_i: C -> C */
  /* typename Problem::I_t */
  auto i = Pb.mix_default();
  
  /* Generating a starting point should be independent in each process,
   * otherwise, we are repeatign the same work in each process!   */
  PRNG prng_elm;

  // ======================================================================== //
  //                                 STEP 2                                   //
  // ------------------------------------------------------------------------ //


  // ======================================================================== //
  //                                 STEP 3                                   //
  // ------------------------------------------------------------------------ //
  bool found_golden_inputs = false;
  
  while (not found_golden_inputs){
    found_golden_inputs = receiver_round<Problem,  Types...>();

    /* Update seeds  */
    /* Clear receiv buffers */
    std::memset(rcv_buf, 0, my_info.nelements_buffer * Pb.C.size);
    std::memset(work_buf, 0, my_info.nelements_buffer * Pb.C.size);
    /* Clear dictionary */
    dict.flush();
    /* repeat! */
  }

  delete[] rcv_buf;
  delete[] work_buf;
}
}


#endif
