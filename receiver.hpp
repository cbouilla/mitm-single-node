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
#include "claw_engine.hpp"
#include "collision_engine.hpp"
#include "mpi_common.hpp"


namespace mitm {

/* Deserialize a message, query the dictionary,
 * 
 */
template <typename Problem, typename DIGEST_T, typename... Types>
bool treat_received_msg(Problem& Pb,
			int const max_nmsgs,
			int const nmsgs,
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

  // Don't change these two without changing them in sender.hpp
  // todo these choices were determined when the dictionary was created!
  using OUT_T = u64; 
  using CHAIN_LENGTH_T = u32;

  // Received buffer:
  // [outputs||chain_lengths||inputs]
  // Index of the first chain_lengths in receive buffer
  size_t const offset_out = 0; 
  size_t const offset_chain_lengths = max_nmsgs * sizeof(OUT_T);
  size_t const offset_inp = offset_chain_lengths + (max_nmsgs * sizeof(OUT_T) );

  // Extracted data from 
  OUT_T digest = 0; 
  CHAIN_LENGTH_T chain_length0 = 0;
  CHAIN_LENGTH_T chain_length1 = 0; 

  // Basic information about the queried triples
  bool matched = false;
  bool found_golden = false;

  for (int msg_i = 0; msg_i < nmsgs; ++msg_i){
    // 1- Extract data
    // 1.a digest
    std::memcpy(&digest,
		&rcv_buf[offset_out + msg_i*sizeof(digest)],
		sizeof(OUT_T));

    // 1.b chain_length
    std::memcpy(&chain_length0,
		&rcv_buf[offset_chain_lengths + msg_i*sizeof(OUT_T)],
		sizeof(CHAIN_LENGTH_T));

    // 1.c input
    Pb.C.deserialize(&rcv_buf[offset_inp + msg_i*Pb.C.size], inp_St);


    // 2- Query the dictionary
    matched = dict.pop_insert(digest,
			      inp_St,
			      chain_length0,
			      *inp1_pt,
			      chain_length1,
			      Pb);

      
    // 3- treat collision if any
    if (matched)
      found_golden = treat_match(Pb,
				 byte_hasher,
				 fn_idx,
				 inp0_pt,
				 out0_pt,
				 chain_length0,
				 inp1_pt, /* todo fix this, todo: what's wrong with it?! */
				 out1_pt,
				 chain_length1,
				 inp_mixed,
				 inp0A, 
				 inp1A,
				 args...); /* for claw args := inp0B, inp1B */
  }
  return found_golden;
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
  // todo find a way to share the counter_size between senders and receivers?
  int counter_size = 2;
  u16 nmsgs = 0; // # triples received 
  MPI_Request request_finish{};
  MPI_Request request{};
  int round_completed = false; // todo check that "false"" maps to zero and "not zero" to 1
  int receive_completed = false;
  bool found_golden = false;
  
  
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
    
    found_golden = treat_received_msg(Pb,
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

    if (found_golden)
      // we found the golden pair, they are inp0A,  inp1A (collision),
      return true; // or  inp0A, inp1B (claw), inp1B is inside args... sorry :( 


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

  

  // ======================================================================== //
  //                                 STEP 1                                   //
  //      Initialization, and  agreeing on mix_function and walk function     //
  // ------------------------------------------------------------------------ //
  u64 mixer_seed;
  u64 byte_hasher_seed;

  /* for MPI tags to block communication when different version of the function 
  * is used. */
  size_t round_number = 0;

  seed_agreement(my_info, mixer_seed, byte_hasher_seed);

  /* Now, mixing function and byte_hasher are the same among all processes */
  PRNG prng_mix{mixer_seed}; /* for mixing function */
  PearsonHash byte_hasher{byte_hasher_seed}; /* todo update PearsonHash to accept a seed */
  /* variable to generate families of functions f_i: C -> C */
  /* typename Problem::I_t */
  typename Problem::I_t fn_idx = Pb.mix_default();
  
  /* Generating a starting point should be independent in each process,
   * otherwise, we are repeatign the same work in each process!   */
  PRNG prng_elm;

  // ======================================================================== //
  //                                 STEP 2                                   //
  //                        receiver, treat, and repeat!                      //
  // ------------------------------------------------------------------------ //


  bool found_golden = false;
  
  while (not found_golden){
    // Main computation step
    found_golden = receiver_round(Pb,
				  my_info,
				  round_number,
				  difficulty,
				  byte_hasher,
				  dict,
				  rcv_buf,
				  fn_idx,
				  inp_St,
				  inp0_pt,
				  inp1_pt,
				  out0_pt,
				  out1_pt,
				  inp_mixed,
				  inp0A,
				  inp1A,
				  args...);
    if (found_golden)
      break; // exit the loop
    
    // Preparation and synchronization
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
