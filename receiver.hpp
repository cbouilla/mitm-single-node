 #ifndef MITM_MPI_RECEIVER
#define MITM_MPI_RECEIVER

#include <cstddef>
#include <cstring>
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

/* Deserialize a message, query the dictionary,
 * 
 */
template <typename Problem, typename DIGEST_T, typename... Types>
bool treat_received_msg(Problem& Pb,
			int nmsgs,
			size_t round_number,
			int const difficulty,
			PearsonHash const& byte_hasher,
			Dict<DIGEST_T, typename Problem::C_t, Problem>& dict,
			u8* rcv_buf,
			typename Problem::I_t& fn_idx,
			typename Problem::C_t& inp_St, // Startign point in chain
			typename Problem::C_t* inp0_pt,
			typename Problem::C_t* inp1_pt,
			typename Problem::C_t* out0_pt,
			typename Problem::C_t* out1_pt,
			typename Problem::C_t& inp_mixed,
			typename Problem::A_t& inp0A,
			typename Problem::A_t& inp1A,
			Types... args)
{

  // why not wrap all of this into a struct
  size_t const offset_inp = 0; // todo correct this
  size_t const offset_out = 0; // todo correct this
  size_t chain_length0 = 0;
  size_t chain_length1 = 0; 
  bool matched = false;
  
  u64 digest = 0; // todo pass it as a parameter?
  for (int msg_i = 0; msg_i < nmsgs; ++msg_i){
    // 1- Deserialize // todo start from here
    Pb.C.deserialize(&rcv_buf[offset_inp + msg_i*Pb.C.size], inp_St);
    
    
    /* todo (critical) change the demos to respect the documentation.
     * Changes:
     * length -> size
     * unserialize -> deserialize
     * treat_collision -> treat_match that calls:
     *  treat_collision (collision_enging) and treat_claw (claw_engine)
     */
    // Get  the digest
    std::memcpy(&digest,
		&rcv_buf[offset_out + msg_i*sizeof(digest)],
		sizeof(digest));
    
    // 2- Query the dictionary
    matched = dict.pop_insert(digest,
			      inp_St,
			      chain_length0,
			      *inp1_pt,
			      chain_length1,
			      Pb);

      
    // 3- treat collision if any
    if (matched)
      continue; // todo continue from here 
    // todo rename treat_collision to treat_match then inside treat_collision or treat_claw
    // repeat for the number of messages received
  }
}


/* One round of computation of receiver  */
template <typename Problem, typename DIGEST_T, typename... Types>
bool receiver_round(Problem& Pb,
		    MITM_MPI_data& my_info,
		    size_t round_number,
		    int const difficulty,
		    PearsonHash const& byte_hasher,
		    Dict<DIGEST_T, typename Problem::C_t, Problem>& dict,
		    u8* rcv_buf,
		    typename Problem::I_t& fn_idx,
		    typename Problem::C_t& inp_St, // Startign point in chain
		    typename Problem::C_t* inp0_pt,
		    typename Problem::C_t* inp1_pt,
		    typename Problem::C_t* out0_pt,
		    typename Problem::C_t* out1_pt,
		    typename Problem::C_t& inp_mixed,
		    typename Problem::A_t& inp0A,
		    typename Problem::A_t& inp1A,
		    Types... args)
{
  // Keep posted to senders when they are done.
  int ndones = 0;
  // todo find a way to make it shared between sender and receiver!
  int counter_size = 2;
  u16 nmsgs = 0; // # triples received 
  MPI_Request request_finish{};
  MPI_Request request{};
  int round_completed = false; // todo check false maps to zero and not zero = 1
  int receive_completed = false;

  
  MPI_Iallreduce(NULL, /* Don't send anything */
		 &ndones, /* receive sum_{done process} (1) */
		 1, // Receive one int
		 MPI_INT,
		 MPI_SUM, // reduction operation
		 my_info.inter_comm, // get the data from senders
		 &request_finish);

  // Listen to received messages.
  while (not round_completed){
    MPI_Irecv(rcv_buf,
	      my_info.msg_size,
	      MPI_UNSIGNED_CHAR,
	      MPI_ANY_SOURCE,
	      ROUND_SND_TAG + round_number,
	      my_info.inter_comm,
	      &request);

    // Get the actual number of messages received
    std::memcpy(&nmsgs,
		&rcv_buf[my_info.msg_size - counter_size],
		counter_size);
    
    treat_received_msg(Pb,
		       nmsgs,
		       round_number,
		       difficulty,
		       byte_hasher,
		       dict,
		       rcv_buf,
		       fn_idx,
		       inp_St, // Startign point in chain
		       inp0_pt,
		       inp1_pt,
		       out0_pt,
		       out1_pt,
		       inp_mixed,
		       inp0A,
		       inp1A,
		       args...);

    /* query those elements in the dictionary, if collision found treat this
    * collision on the spot (call a function that calls a function to treat
    * a collision since we may leave collision treatmenets to another processes)
    */

    /* Assume we only check if all senders completed, if not, then check if 
     * a receive was completed. This lead to a deadlock!
     * Example: Assume all senders completed, except one that is going to
     *          send to another receiver. If we wait for MPI_Irecv to complete
     *          we would wait forever since all senders are completed except 
     *          one not assigned to this receiver.
     */
    while (not (round_completed or receive_completed) ){
      MPI_Test(&request, &receive_completed, MPI_STATUS_IGNORE);
      MPI_Test(&request_finish, &round_completed, MPI_STATUS_IGNORE);
      if (request_finish)
	return false; // we have not found golden collision yet!
    }
  }
  // If that happens to be everyone, then exit the function.

}



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
