#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER

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
/******************************************************************************/

// -------------------------------------------------------------------------- //
// ---------------------------------- SENDER --------------------------------- //

/* Create two type:
 * 1) working buffers where we book t bytes for each receiver.
 * 2) sending buffers that are ready to be sent to specific receivers.
 * Note: OUT_T is the digest of the output, we limit ourselves to 64bit digest.
 *       CHAIN_LENGTH_T choosen to be u32 since we expect difficulty to be less 32bits
 */ // inputs don't need to be templated since they are sequence of bytes.
//template <typename OUT_T = u64, typename CHAIN_LENGTH_T = u32>
struct sender_buffers {

  using OUT_T = u64;
  using CHAIN_LENGTH_T = u32;
  // Assume there are `n` receivers, per receiver we send `m` elements, each of
  // size `l`.
  // Then to get the ith element that will be sent to receiver `r`:
  u8* inputs;   // ith input  <-  inputs[(r * m * l) +  i * l]
  u64* outputs; // ith output <- outputs[(r * m)     +  i]
  u32* chain_lengths; // ith chain length <- chain_length[(r * m)  +  i]

  // ith element <- #elements stored in reciever `r` buffer. 
  // We use u16 since nelements to be sent  are ~1k at a time < 2^16
  u16* counters; 
  size_t const counter_size = sizeof(u16);
  
  // When a receiver `r' buffers (inputs, outputs, chain_lengths) gets filled,
  // it will be copied to `snd` buffer. Then, MPI will upload `snd` to the
  // receiver `r`.
  // In `snd` all values are chained as the following
  //  snd =  outputs || chain_length || inputs || nelements,
  // i.e.
  // ith output <-  cast_u64(snd[8*i]), recall we send digests as the outputs
  // ith chain_length <- cast_u32( snd[nelements*8 + i*4] )
  // ith input <- snd[nelements*8 + nelements*4 + i*m]
  u8* snd; 
  // Maximum number of elements in a single send.
  size_t const nelements_snd_max;
  /* Regardless of the number of elements in the buffer, always send fixed  */
  size_t const snd_size; // amount of bytes.
  size_t const nreceivers; 
  size_t const inp_size;
  size_t const inp_out_chain_size;
  // How many bytes for counter that counts how many triples in the send buffer

  
  /* sender has: 1) working buffers. 2) sending buffers ready to be sent  */
  sender_buffers(size_t const nelements,
		 size_t const nreceivers,
		 size_t const inp_size)
    // The last 4 bytes for #triples in the buffer.
    : nelements_snd_max(nelements),
      snd_size(nelements * (sizeof(OUT_T) + sizeof(u32) + inp_size) + counter_size),
      nreceivers(nreceivers),
      inp_size(inp_size),
      inp_out_chain_size(sizeof(OUT_T) + sizeof(u32) + inp_size)
  {
    inputs  = new u8[nreceivers * nelements_snd_max * inp_size];

    // on the other hand, we are only going send the digests.
    //outputs = new u64[nreceivers * nelements ];
    outputs = new OUT_T[nreceivers * nelements_snd_max];

    // How many times we applied f(inp) to get the corresponding output
    // chain_lengths = new u32[nreceivers * nelements ];
    chain_lengths = new CHAIN_LENGTH_T[nreceivers * nelements_snd_max ];
 
    // snd = new u8[nelements*sizeof(u64)  // outputs of one receiver 
    // 		 + nelements*sizeof(u32) // chains lengths of 1 receiver
    // 		 + nelements*element_size]; // inputs of receiver
    snd = new u8[snd_size]; // + 4 for ntriples in buffer
    

    counters = new u16[nreceivers];

    // Not sure, but I feel it's better to map virtual memory to ram right away!
    inputs[0]        = 0;
    outputs[0]       = 0;
    chain_lengths[0] = 0;
    snd[0]           = 0;
    
    // This is the only buffer that requires to be zero. 
    std::memset(counters, 0, nreceivers*sizeof(u16));
  }

  ~sender_buffers() /* free up the memory */
  {
    delete[] inputs ;
    delete[] outputs;
    delete[] chain_lengths;
    delete[] snd;
  }

  /* Clear buffers for next use */
  void clear()
  { 
    std::memset(counters, 0, nreceivers*sizeof(u16));  
  }

  /* Reset buffer i counter to zero, and clear and other data if necessary */
  void clear_buf(int id)
  {
    counters[id] = 0;
  }

  /* Put all data found in receivers number i buffers (inputs, outputs, chains)
   * into snd_buffer. Use counters[i] to see how many elements to put in the 
   */
  void copy_receiver_i_to_snd(int id)
  {
    // Recall: snd = [outputs] || [chain lengths] || [inputs]
    // where each array has fixed size |element| * nelements_max_snd regardless
    // of the number of actual elements.
    
    // snd_size = static_cast<u64>(counters[id]) * inp_out_chain_size
    //          + counter_size;

    // todo: check this sizeof(u64) * counters[id] is on u64;
    //    size_t offset = sizeof(OUT_T) * counters[id];  // wrong
    
    //size_t offset = 0; 
    std::memcpy(&snd[0],
		outputs,
		sizeof(OUT_T) * counters[id]);

    // It is always constant regardless of #elements stored.
    size_t const offset_chain = sizeof(OUT_T) * nelements_snd_max;
    std::memcpy(&snd[offset_chain],
		chain_lengths,
		sizeof(CHAIN_LENGTH_T) * counters[id]);

    size_t const offset_out = (sizeof(CHAIN_LENGTH_T) * nelements_snd_max)
                            + offset_chain;
    
    std::memcpy(&snd[offset_chain],
		inputs,
		inp_size * counters[id]);

    // Finally add how many elements were copied.
    size_t const offset_ctr = offset_out + offset_chain + (inp_size * counters[id]);
    std::memcpy(&snd[offset_ctr],
		&counters[id],
		counter_size);
  }

  int to_which_receiver(u64 digest)
  {
    return 0; // todo fix this!
  }
};


/* Fill receivers buffers until one of them gets filled.
 * RETURN: receiver id whose buffer is full.
 */
template<typename Problem,  typename... Types>
int fill_buffers(Problem& Pb,
		 int const difficulty,
		 PRNG& prng_elm,
		 PearsonHash const& byte_hasher,
		 sender_buffers& buffers, // constant pointer
		 typename Problem::I_t& func_idx,
		 typename Problem::C_t& inp_St, // Startign point in chain
		 typename Problem::C_t* inp0_pt,
		 typename Problem::C_t* inp1_pt,
		 typename Problem::C_t* out0_pt,
		 typename Problem::C_t* out1_pt,
		 typename Problem::C_t& inp_mixed,
		 typename Problem::A_t& inp0A,
		 typename Problem::A_t& inp1A,
		 Types... args)/* two extra arguments if claw problem */
{
  size_t chain_length = 0;
  int found_dist = -1; // Was the search for a disitinguished point successful?
  bool is_full_buffer = false;
  size_t out_digest = 0;
  size_t receiver_id = 0;
  u16 buff_nelements = 0;
  
  while(not is_full_buffer){
    Pb.C.randomize(inp_St, prng_elm); // Get new fresh element.
    Pb.C.copy(inp_St, *inp0_pt); // Copy the starting point to buffer that changes

    found_dist = generate_dist_point(Pb,
				     byte_hasher,
				     func_idx, /* permutation number of f's input */
				     chain_length,
				     difficulty,
				     inp0_pt,
				     out0_pt,
				     inp_mixed,
				     inp0A,
				     args...); /* for claw args := inp0B, inp1B */
    if (found_dist){
      out_digest = Pb.C.hash(*out0_pt);
      // to which receiver?
      receiver_id = (out_digest >> difficulty) % buffers.nreceivers;
      // increment nelements counter
      buff_nelements = buffers.counters[receiver_id]; 
      // Copy the input, output, chain_length to receiver_id buffers
      buffers.outputs[buff_nelements] = out_digest;
      buffers.chain_lengths[buff_nelements] = chain_length;
      Pb.C.serialize(inp_St, // from here to the location below // todo use Pb.C.size
		     &buffers.inputs[buff_nelements * buffers.inp_size]
		     );
      // record the new number of elements in receiver_id buffer 
      ++buffers.counters[receiver_id]; 

      if (buffers.counters[receiver_id] >= buffers.nelements_snd_max)
	is_full_buffer = true; // break from the loop
    }
  }

  return receiver_id; // last one whose buffer became full.
}

/* Send x messages to receivers, after that update, return was golden collision found?  */
template<typename Problem,  typename... Types>
bool sender_round(Problem& Pb,
		  MITM_MPI_data& my_info,
		  size_t round_number, // used to differetiate between different function updates
		  int const difficulty,
		  PearsonHash const& byte_hasher,
  		  sender_buffers& buffers, // constant pointer
		  typename Problem::I_t& func_idx,
		  typename Problem::C_t& inp_St, // Startign point in chain
		  typename Problem::C_t* inp0_pt,
		  typename Problem::C_t* inp1_pt,
		  typename Problem::C_t* out0_pt,
		  typename Problem::C_t* out1_pt,
		  typename Problem::C_t& inp_mixed,
		  typename Problem::A_t& inp0A,
		  typename Problem::A_t& inp1A,
		  Types... args)/* two extra arguments if claw problem */

{
  MPI_Request request;

  u64 elm_seed = read_urandom<u64>();
  /* Generating a starting point should be independent in each process,
   * otherwise, we are repeatign the same work in each process!   */
  PRNG prng_elm{elm_seed};

  
  int buf_id = -1;
  size_t beta = 10; // Wiener&Oorschost say 
  size_t nsends = beta * (my_info.nelements_mem_receiver * my_info.nsenders);
  
  for (size_t i = 0; i < nsends; ++i){ // tood x is the number of times to fill buffer

    /* Fill buffers until one becomes full */
    buf_id =  fill_buffers(Pb,
			   difficulty,
			   prng_elm,
			   byte_hasher,
			   func_idx,
			   buffers,
			   inp_St, inp0_pt, inp1_pt, out0_pt, out1_pt,
			   inp_mixed, inp0A, inp1A, args...);
    buffers.clear_buf(buf_id); // reset the buffer counter.

    if (i > 0) // skip MPI_Wait in the first round
      MPI_Wait(&request, MPI_STATUSES_IGNORE); // wait for the previous send.

    // todo copy buffers to snd buffer!!!
    buffers.copy_receiver_i_to_snd(buf_id);
    
    /* Send the buffer that is already full  */
    MPI_Isend(buffers.snd,
	      buffers.snd_size,
	      MPI_UNSIGNED_CHAR,
	      buf_id,
	      ROUND_SND_TAG + round_number, // todo don't forget to add round number to the tag. Also, check ROUND_SND_TAG is the largest tag. // todo edit common_mpi.hpp
	      my_info.inter_comm,
	      &request);
  }

  /* send the remaining elements in the buffer */
  for (size_t id = 0; id < buffers.nreceivers; ++id){
    buffers.copy_receiver_i_to_snd(id);
    
    /* Send the buffer that is already full  */
    MPI_Isend(buffers.snd,
	      buffers.snd_size,
	      MPI_UNSIGNED_CHAR,
	      id,
	      ROUND_SND_TAG + round_number, // todo don't forget to add round number to the tag. Also, check ROUND_SND_TAG is the largest tag. // todo edit common_mpi.hpp
	      my_info.inter_comm,
	      &request);
  }

  // todo how sender would know a golden collision was found?
  return false; /* no golden collision was found */
}


template<typename Problem,  typename... Types>
void sender(Problem& Pb,
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
  /* Pseudo-random number generators for elements */


  
  /* ----------------------------- sender buffers ---------------------------- */
  sender_buffers buffers(my_info.nelements_buffer,
			 my_info.nreceivers,
			 Pb.C.size);
  
  

  
  
  // ======================================================================== +
  //                                 STEP 1                                   +
  // ------------------------------------------------------------------------ +
  //      Initialization, and  agreeing on mix_function and walk function     +
  // ------------------------------------------------------------------------ +
  /* Top seeds to be agreed among all processes */
  u64 seed1;
  u64 seed2;

  seed_agreement(my_info, seed1, seed2);
  // ------------------------------------------------------------------------+
  // Create my local PRNG that agrees with everyone globally.
  /* Then we generate the next seeds using prng (deterministic process) */
  PRNG mixer_seed_gen{seed1};
  u64 mixer_seed = mixer_seed_gen.rand();

  PRNG byte_hasher_seed_gen{seed2};
  u64 byte_hasher_seed = byte_hasher_seed_gen.rand();

  /* Now, mixing function and byte_hasher are the same among all processes */
  PRNG prng_mix{mixer_seed}; /* for mixing function */
  PearsonHash byte_hasher{byte_hasher_seed};
  
  /* variable to generate families of functions f_i: C -> C */
  /* typename Problem::I_t */
  auto i = Pb.mix_default();
  // -------------------------------------------------------------------------+

  PRNG prng_elm;

  // ======================================================================== +
  //                                 STEP 2                                   +
  // ------------------------------------------------------------------------ +

   while (true){
     sender_round(Pb,
		  my_info,
		  difficulty,
		  byte_hasher,
		  buffers,
		  inp_St,
		  inp0_pt, inp1_pt, out0_pt, out1_pt, inp_mixed,
		  inp0A, inp1A,
		  args...);
     

     /* Update seeds  */
     mixer_seed = mixer_seed_gen.rand();
     byte_hasher_seed = byte_hasher_seed_gen.rand();

     /* update mix_function */
     i = Pb.mix_sample(prng_mix);
     
     /* update byte_hasher (claw) */
     byte_hasher.update_table(byte_hasher_seed);
     
     /* Clear buffers */
     buffers.clear();
     
    /* repeat! and good luck! */
  }
}

}
#endif
