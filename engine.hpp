#ifndef MITM_ENGINE
#define MITM_ENGINE

/*****************************************************************************/
// GLOBAL VARIABLES FOR DEBUGGING THE NUMBER OF COLLISIONS BEFORE THE GOLDEN
// 0xdeadbeef tag for debugging 
/* we saw the golden input, and in the chain C -f/g-> C it uses f: A -> C */
#ifdef CLAW_DEBUG
bool found_golden_A_and_use_f = false;
/* we saw the golden input, and in the chain C -f/g-> C it uses g: B -> C */
bool found_golden_B_and_use_g = false;
/**/
#endif

#ifdef COLLISION_DEBUG
bool found_1st_golden_inp = false;
bool found_2nd_golden_inp = false;
#endif 
/*****************************************************************************/
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




/// When a function, func, behavior's depends if the problem is a claw or a
///  collision, we use veriadic template to have two different implementation:
/// template<_, typename... Types>  /* A veriadice template */
/// func(_, Types... args)
/// Then in `claw_engine.hpp` would have a template specialization of func
/// with more arguments than the specialization of func found in
/// `collisions_engine.hpp`
/// List of functions that are different for claw and collision problem:
/// 1) `iterate_once` claw_engine.hpp introduces 2 more arguments.

namespace mitm {

/******************************************************************************/
/* The two functions below are implemented in `claw_engine.hpp` and
 *  `collision_engine.hpp`
 */
  
/*
 * Defines how to walk in a sequence. Use functions overloading to differentiate
 * between a claw's iteration and a collisions's iteration. Specfically:
 * When it's a claw problem, then it would have 2 extra arguments, namely:
 * 1)  inp0B, and inp1B.
 * When it's a collision:
 * args is empty
 * 
 */
// template <typename Problem, typename... Types>
// void iterate_once(Problem &Pb,
// 		  typename Problem::I_t& i, /*permutation number of f's input*/
// 		  typename Problem::C_t& inp,
// 		  typename Problem::C_t& out,
// 		  typename Problem::C_t& inp_mixed,
// 		  typename Problem::A_t& inpA, /* scratch buffer */
// 		  Types... args) /* for claw args := inp0B, inp1B */
// {
//   /* C++ always prefers more specialized templates. 
//    * For the claw code, see `iterate_once` in `claw_engine.hpp`
//    * For the collsions code, see `iterate_once` in `collision_engine.hpp`
//    */
//   std::cerr << "Internal Error: should not use general implementation of `iterate_once`!\n";
//   std::terminate(); /* Never use this implementation! */
// }



/*
 * Given two collisions in the modified function 'f: C_t -> C_t send them to
 * their original domain: A_t or B_t (if it's a claw), and test if they make
 * THE golden_collision. The implementation of this function is found in
 * `claw_engine.hpp` or `collision_engine.hpp`
 * In case of claw: an additional arguments:
 * 1) inpB_pt
 */
// template <typename Problem, typename... Types >
// bool treat_collision(Problem& Pb,
// 		     typename Problem::I_t& i,
// 		     typename Problem::C_t*& inp0_pt,
// 		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
// 		     const u64 inp0_chain_len,
// 		     typename Problem::C_t*& inp1_pt,
// 		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
// 		     const u64 inp1_chain_len,
// 		     typename Problem::C_t& inp_mixed,
// 		     typename Problem::A_t& inp0_A,
// 		     typename Problem::A_t& inp1_A,
// 		     Types... args) /* for claw args := inp0B, inp1B */
// {
//   std::cerr << "Internal Error: should not use general implementation of if `treate_collision`!\n";
//   std::terminate(); /* Never use this implementation! */
//}

/******************************************************************************/
/*                                     MPI                                    */


enum process_type{ SENDER,  RECEIVER };
enum MITM_MPI_TAGS {INTERCOM_TAG, ROUND_SND_TAG}; // to be extended ...


struct MITM_MPI_data{
  MPI_Comm global_comm; // Global communicator. Kept just in case, todo to be removed!
  MPI_Comm local_comm;  // Processes that do the same thing, e.g. senders. 
  MPI_Comm inter_comm;  // splitted into two: senders and receivers.
  size_t const nelements_buffer = 1000; // #elements stored in buffer during send/receive
  size_t nelements_in_mem; // How many elements can be stored across all dictionaries.
  int nprocesses;      // Total number of processes.
  int nsenders;         // Total number of senders .
  int nreceivers;       // Total number of receivers.
  int nnodes; // How many nodes are there? 
  int nsenders_node;    //  nsenders on my node (counts me if i'am a sender).
  int nreceivers_node;  //  nreceivers on my node (counts me if i'am a receiver).
  int my_rank_intercomm; // My rank among the group of the same color
  int my_rank_global;
  process_type my_role; // am I a sender or a receiver.
};

    
/* Get the number of number of physical nodes that runs the program.
 * It requires MPI 3 or above.
 * source: https://stackoverflow.com/a/34118174 */
int get_nodes_count(MPI_Comm comm)
{
  int rank, is_rank0, nodes;
  MPI_Comm shmcomm;

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
		      MPI_INFO_NULL, &shmcomm);
  
  MPI_Comm_rank(shmcomm, &rank);
  is_rank0 = (rank == 0) ? 1 : 0;
  MPI_Allreduce(&is_rank0, &nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Comm_free(&shmcomm);

  return nodes;  
}


/* Return how many bytes of RAM can a receiver allocate,
 * if the available space is shared. */
size_t get_nbytes_per_receivers(MITM_MPI_data& my_info)
{
  return 0;
}

/* Initializes mpi processes for mitm,  and splits senders and receivers into
 * seperate communicators. Aslo, gives the number of total numbers of senders
 *  and receivers, `nsenders` and `nreceivers` respectively,  and the number
 * senders and receivers that are on the same node, `nsenders_node` & 
 * `nreceivers_node` respectively. This allows us to use nodes with different
 * specs without affecting the code.
 * 
 */
MITM_MPI_data MITM_MPI_Init(int nreceivers, int argc=0, char** argv = NULL)
{
  MITM_MPI_data my_info{};
  

  /* Basic MPI Boilerplate */
  MPI_Init(&argc, &argv); /* Provided argc, and argv to silence C++ errors */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_info.my_rank_global);
  MPI_Comm_size(MPI_COMM_WORLD, &my_info.nprocesses);

  /* Two disjoint groups: receivers, and senders */
  int is_receiver = (my_info.my_rank_global < nreceivers);

  /* Create an inter-comm that seperates senders and receivers */
  MPI_Comm_split(MPI_COMM_WORLD, /* Communicator to be splitted */
		 is_receiver, /* create comm with processes that share this value  */
		 my_info.my_rank_global, 
		 &my_info.inter_comm /* New communicator. */
		 ); /* sad face for the seperation between processes */


  /* Create communicator between the two clans of processes. */
  if (is_receiver) {
    /*  Receivers establish their communication channel with senders */
    MPI_Intercomm_create(my_info.local_comm, /* my tribe */
			 0, /* Local name of the sheikh of my tribe */
			 MPI_COMM_WORLD, /* Where to find the remote leader. */
			 nreceivers, /* Global Name/rank of senders leader. */
			 INTERCOM_TAG, /* This tag is unique to the creation of intercomm */
			 &my_info.inter_comm); /* <- put down these info here */
    my_info.my_role = RECEIVER;
    
  } else {
    /* senders establish their communication channel with receivers */
    MPI_Intercomm_create(my_info.local_comm, /* my tribe */
			 0, /* Local name of the leader of my tribe */
			 MPI_COMM_WORLD, /* Where to find the remote leade. */
			 0, /* Global name/rank of receivers leader. */
			 INTERCOM_TAG, /* This tag is unique to the creation of intercomm */
			 &my_info.inter_comm); /* <- put down these info here */

    my_info.my_role = SENDER;
  }

  /* Fill the rest of the details  */
  my_info.nsenders = my_info.nprocesses - nreceivers;
  my_info.nreceivers = nreceivers;
  my_info.nnodes = get_nodes_count(MPI_COMM_WORLD);
  my_info.nsenders_node = my_info.nsenders / my_info.nnodes;
  my_info.nreceivers_node = my_info.nreceivers / my_info.nnodes;
  
  MPI_Comm_rank(my_info.local_comm, &my_info.my_rank_intercomm);
  

  /* We can say how many elements can be stored in memory */
  // todo formulas are incorrect!
  size_t nbytes_avaiable = get_available_memory();
  size_t sender_needed_bytes = 0;
  size_t dict_one_element_size = 1;
  // memory available to all receivers per node 
  size_t nelements_ram_receiver = (nbytes_avaiable - sender_needed_bytes);
  // nelements that will be stored in all receivers per node
  nelements_ram_receiver = nelements_ram_receiver / dict_one_element_size;
  // nelements per receiver
  nelements_ram_receiver = nelements_ram_receiver / my_info.nreceivers_node;

  my_info.nelements_in_mem = nelements_ram_receiver;

  return my_info;
}


/******************************************************************************/
// Model of send receive:
// 1) Create sender and receivers. Also, determine the number of elements
//    the receivers have.
// 2) Agree on seed by brodacasting the seed from rank 0 to everyone
// 3) A sender fill receivers buffer and send as soon as a buffer is complete.
//    Also, send the number of elements along with the message. This helps to
//    say if sender has finished its round (all senders would do the same).
//    The tag used is snd_tag + i where i is the round number to block
//    from processing elements from the next round.
// 4) recievers receive messages. increase a counter if the number of received
//    is less than the maximum. update seeds when the counter is equal the
//    equal the number of senders
// 
// 


/******************************************************************************/

// -------------------------------------------------------------------------- //
// ---------------------------------- SENDER --------------------------------- //

/* Create two type:
 * 1) working buffers where we book t bytes for each receiver.
 * 2) sending buffers that are ready to be sent to specific receivers. 
 */
struct sender_buffers {
  // we could have an array of triple (inp, out, digest).
  // Assume there are `n` receivers, per receiver we have `m` elements, each of
  // size `l`. Then to get the ith element that will be sent to receiver `k`:
  u8* inputs;   // ith input  <-  inputs[(k * m * l) +  i * l]
  u64* outputs; // ith output <- outputs[(k * m)     +  i]
  u32* chain_lengths; // ith chain length <- chain_length[(k * m)  +  i]

  // ith element <- #elements stored in reciever i buffer. We use u16 since
  // we are going to send ~1k at a time, i.e. nelements < 2^10 
  u16* counters; 
  
  // One receiver buffers to be sent where all values chained as 
  //  snd =  outputs || chain_length || inputs, i.e.
  // ith output <-  cast_u64(snd[8*i]),
  // ith chain_lenght <- cast_u32( snd[nelements*8 + i*4] )
  // ith input <- snd[nelements*8 + nelements*4 + i*m]
  u8* snd; 

  size_t const snd_size_max;
  size_t const element_size; /* we store this information twice! */
  /* Au cas où, get a segmentation error if it used incorrectly. */
  size_t snd_size = -1;
  size_t counters_size;

  /* sender has: 1) working buffers. 2) sending buffers ready to be sent  */
  sender_buffers(size_t nelements, size_t nreceivers, size_t element_size)
    : snd_size_max(nelements * (sizeof(u64) + sizeof(u32) + element_size)),
      element_size(element_size),
      counters_size(nreceivers)
  {

    inputs  = new u8[nreceivers * nelements * element_size];

    // on the other hand, we are only going send the digests.
    outputs = new u64[nreceivers * nelements ];

    // How many times we applied f(inp) to get the corresponding output
    chain_lengths = new u32[nreceivers * nelements ];
 
    // snd = new u8[nelements*sizeof(u64)  // outputs of one receiver 
    // 		 + nelements*sizeof(u32) // chains lengths of 1 receiver
    // 		 + nelements*element_size]; // inputs of receiver
    snd = new u8[snd_size_max + 4]; // + 4 for ntriples in buffer 
    

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
    std::memset(counters, 0, counters_size*sizeof(u16));  
  }

  /* Put all data found in receivers number i buffers (inputs, outputs, chains)
   * into snd_buffer. Use counters[i] to see how many elements to put in the 
   */
  void copy_receiver_i_to_snd(int id, int msg_size)
  {
    snd_size = static_cast<u64>(counters[id])
             * (sizeof(u64) + sizeof(u32) + element_size);

    // todo: check this sizeof(u64) * counters[id] is on u64;
    size_t offset = sizeof(u64) * counters[id]; 
    std::memcpy(&snd[0], outputs, offset);
    std::memcpy(&snd[offset], chain_lengths, sizeof(u32) * counters[id]);

    offset += sizeof(u32) * counters[id];
    std::memcpy(&snd[offset], inputs, element_size * counters[id]);

    // todo add msg_size at the end of the buffer
    
  }
};

    
template<typename Problem,  typename... Types>
int sender_fill_buff(Problem& Pb,
		     int const difficulty,
		     PearsonHash const& byte_hasher,
		     sender_buffers& buffers, // constant pointer 
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
  
}

/* Send x messages to receivers, after that update, return was golden collision found?  */
template<typename Problem,  typename... Types>
bool sender_round(Problem& Pb,
		  MITM_MPI_data& my_info,
		  int const difficulty,
		  PearsonHash const& byte_hasher,
  		  sender_buffers& buffers, // constant pointer 
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

  
  int buf_id =  sender_fill_buff(Pb,
				 difficulty,
				 byte_hasher,
				 buffers,
				 inp_St, inp0_pt, inp1_pt, out0_pt, out1_pt,
				 inp_mixed, inp0A, inp1A, args...);

  size_t beta = 10; // Wiener&Oorschost say 
  size_t nsends = beta * (my_info.nelements_in_mem / my_info.nsenders);
  
  for (size_t i = 0; i < nsends; ++i){ // tood x is the number of times to fill buffer
    /* Send the buffer that is already full  */
    MPI_Isend(buffers.snd,
	      buffers.snd_size,
	      MPI_UNSIGNED_CHAR,
	      buf_id,
	      ROUND_SND_TAG,
	      my_info.inter_comm,
	      &request);

    /* Fill buffers until one becomes full */
    buf_id =  sender_fill_buff(Pb,
			       difficulty,
			       byte_hasher,
			       buffers,
			       inp_St, inp0_pt, inp1_pt, out0_pt, out1_pt,
			       inp_mixed, inp0A, inp1A, args...);

    /* Check that the previous sending was completed */
    MPI_Wait(&request, MPI_STATUSES_IGNORE);
  }

  /* send the remaining elements in the buffer */
  // todo complete this section
  
  return false; /* no golden collision was found */
}

/* One round of computation of receiver  */
template <typename Problem, typename... Types> bool receiver_round() {}


/* Everyone should agree on the same seed at the beginning.
 * Update the arguments: mixer_seed, byte_hasher_seed */
void seed_agreement(MITM_MPI_data& my_info,
		    u64& mixer_seed,
		    u64& byte_hasher_seed)
{
  u64 seeds[2] = {0};
  int emitter_rank = my_info.nreceivers;

  // if my intercomm rank = 0, i.e. I am the leader
  // I am resposible for sending the initial seed.
  if (my_info.my_rank_global == emitter_rank){
     seeds[0] = read_urandom<u64>();
    seeds[1] = read_urandom<u64>();
  }
  
  /* All sender should agree on mixing function in each round */  
  /* Broadcast seeds to everyone */
  MPI_Bcast(seeds, 2, MPI_UINT64_T, emitter_rank, my_info.global_comm);

  /* Write the received seeds to my memory */
  mixer_seed = seeds[0];
  byte_hasher_seed = seeds[1];

  /* It was a nice exchange, thank you! */
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
  
/******************************************************************************/


inline bool is_distinguished_point(u64 digest, u64 mask)
{  return (0 == (mask & digest) ); }

template<typename C_t>
inline void swap_pointers(C_t*& pt1,
                          C_t*& pt2){
    /// pt1 will point to what pt2 was pointing at, and vice versa.
    C_t* tmp_pt = pt1;
    pt1 = pt2;
    pt2 = tmp_pt;
}


/*
 * Given an input, iterate functions either F or G until a distinguished point
 * is found, save the distinguished point in out_pt and output_bytes
 * Then return `true`. If the iterations limit is passed, returns `false`.
 */
template<typename Problem, typename... Types>
bool generate_dist_point(Problem& Pb,
			 PearsonHash& byte_hasher,
			 typename Problem::I_t& i,
			 u64& chain_length, /* write chain lenght here */
			 const i64 difficulty, /* #bits are zero at the beginning */
			 typename Problem::C_t*& inp_pt,
			 typename Problem::C_t*& out_pt,
			 typename Problem::C_t& inp_mixed,
			 typename Problem::A_t& inpA,
			 Types... args) /* for claw args := inp0B, inp1B */
{
  const u64 mask = (1LL<<difficulty) - 1;
  u64 digest = 0;
  bool found_distinguished = false;
  
  /* The probability, p, of NOT finding a distinguished point after the loop is
   * Let: theta := 2^-d
   * difficulty, N = k*2^difficulty then,
   * p = (1 - theta)^N =>  let ln(p) <= -k
   */
  constexpr u64 k = 40;
  for (u64 j = 0; j < k*(1LL<<difficulty); ++j){
    /* uses claw's iterate_once if args... is not empty, otherwise collisions'*/
    /* for claw args := inp0B, inp1B */
    iterate_once(Pb, byte_hasher,i, *inp_pt, *out_pt, inp_mixed, inpA, args...); 
    ++chain_length;

    /* we may get a dist point here */
    digest = Pb.C.hash(*out_pt);
    found_distinguished = is_distinguished_point(digest, mask);

    /* no need to continue if a distinguished is found  */
    if (found_distinguished){ return true; }

    /* swap inp and out */
    swap_pointers(inp_pt, out_pt);
  }

  return false; /* no distinguished point were found */
}



  
/*
 * A sanity check that verifies the user has provided serialization
 * with its inverse.
 */
template <typename Problem>
bool is_serialize_inverse_of_unserialize(Problem Pb, PRNG& prng)
{
  /// Test that unserialize(serialize(r)) == r for a randomly chosen r
  using C_t = typename Problem::C_t;
  const size_t length = Pb.C.length;

  C_t orig{};
  C_t copy{};
  u8 serial[length];
  size_t n_elements = Pb.C.n_elements;
  const size_t n_tests = std::min(n_elements, static_cast<size_t>(1024));

  for(size_t i = 0; i < n_tests; ++i){
    /* */
    Pb.C.randomize(orig, prng);
    Pb.C.serialize(orig, serial);
    Pb.C.unserialize(serial, copy);

    if (not Pb.C.is_equal(copy, orig)){
      std::cerr << "Error at testing unserial(serial(x)) == x \n";
      std::terminate();
    }
  }
  return true;
}


/* Given two inputs that lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * add a drawing to illustrate this.
 */
template<typename Problem, typename... Types>
bool walk(Problem& Pb,
	  PearsonHash& byte_hasher,
	  typename Problem::I_t& i,
	  u64 inp0_chain_len,
	  typename Problem::C_t*& inp0_pt,
	  typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
	  u64 inp1_chain_len,
          typename Problem::C_t*& inp1_pt,
	  typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
	  typename Problem::C_t& inp_mixed,
	  typename Problem::A_t& inpA,
	  Types... args)
{
  /****************************************************************************+
   *            walk the longest sequence until they are equal                 |
   * Two chains that leads to the same distinguished point but not necessarily |
   * have the same length. e.g.                                                |
   *                                                                           |
   * chain1: ----------------x-------o                                         |
   *                        /                                                  |
   *          chain2: ------                                                   |
   *                                                                           |
   * o: is a distinguished point                                               |
   * x: the collision we're looking for                                        |
   *                                                                           |   
   ****************************************************************************/

  /* Both sequences need at least `len` steps to reach disitinguish point. */
  size_t const len = std::min(inp0_chain_len, inp1_chain_len);

  /* move the longest sequence until the remaining number of steps is equal */
  /* to the shortest sequence. */
  for (; inp0_chain_len > inp1_chain_len; --inp0_chain_len){
    /* for claw args := inp0B, inp1B */
    iterate_once(Pb, byte_hasher, i, *inp0_pt, *out0_pt, inp_mixed, inpA, args...);
    swap_pointers(inp0_pt, out0_pt);
  }
  
  for (; inp0_chain_len < inp1_chain_len; --inp1_chain_len){
    /* for claw args := inp0B, inp1B */
    iterate_once(Pb, byte_hasher, i, *inp1_pt, *out1_pt, inp_mixed, inpA, args...);
    swap_pointers(inp1_pt, out1_pt);
  }


  
  /*****************************************************************************/
  /* now both inputs have equal amount of steps to reach a distinguished point */
  /* both sequences needs exactly `len` steps to reach distinguished point.    */
  for (size_t j = 0; j < len; ++j){
    /* walk them together and check each time if their output are equal     */
    /* return as soon equality is found. The equality could be a robinhood. */

    /* get the outputs of the curren inputs (for claw args... := inp0B, inp1B) */
    iterate_once(Pb, byte_hasher, i, *inp0_pt, *out0_pt, inp_mixed, inpA, args...);
    iterate_once(Pb, byte_hasher, i, *inp1_pt, *out1_pt, inp_mixed, inpA, args...);

    /* First, do the outputs collide? If yes, return true and exit. */
    if(Pb.C.is_equal( *out0_pt, *out1_pt )){
      /* inp0 & inp1 contain  input before mixing, we need to fix this. Maybe yes, maybe no, keep an eye on this comment  */
      return true; /* They are equal */
    }

    /* Move the inputs one step further. The next input is the current output,
     * thus let inp_pt points to the current input data. we don't care what are
     * the data out_pt points to since it will be overwritten by `iterate_once`.
     */
    swap_pointers(inp0_pt, out0_pt);
    swap_pointers(inp1_pt, out1_pt);

  }
  return false; /* we did not find a common point */
}


// template <typename Problem>
// void print_collision_information(typename Problem::C_t& inp0,
// 				 typename Problem::C_t& inp1,
// 				 typename Problem::C_t& out0,
// 				 typename Problem::C_t& out1,
// 				 Problem& Pb)
// {
//   bool real_collision = Pb.C.is_equal(out0, out1);
//   std::cout << "\n++++++++++++++++++++++++++++++++++++++++\n"
// 	    << "Found golden Pair !\n"
// 	    << "inp0 = " << inp0 << "\n"
// 	    << "out0 = " << out0 << "\n"
// 	    << "inp1 = " << inp1 << "\n"
// 	    << "out1 = " << out1 << "\n"
// 	    << "out0 == out1? " << real_collision  << "\n"
// 	    << "++++++++++++++++++++++++++++++++++++++++\n";

// }


/* return a string that says "claw" or "collision" based on the number of
 * arguments */
template <typename... Types>
std::string is_claw_or_collision_problem(Types... args)
{
  if (sizeof...(args) == 2)
    return std::string("collision");

  if (sizeof...(args) == 4)
    return std::string("claw");
  
}

/* The basic search finds collision(s) of
 * iterate_once: C -> C, these collisions can be converted into a solution
 * to the problem we are treating, claw or a collision.
 *  
 */
template<typename Problem,  typename... Types>
void search_generic(Problem& Pb,
		    size_t nbytes_memory,
		    int difficulty,
		    typename Problem::C_t& inp_St, // Startign point in chain
		    typename Problem::C_t* inp0_pt,
		    typename Problem::C_t* inp1_pt,
		    typename Problem::C_t* out0_pt,
		    typename Problem::C_t* out1_pt,
		    typename Problem::C_t& inp_mixed,
		    typename Problem::A_t& inp0A,
		    typename Problem::A_t& inp1A,
		    Types... args) /* two extra arguments if claw problem
				    * 1) inp0B passed as a reference
				    * 2) inp1B passed as a reference */
{
  /* Pseudo-random number generators for elements */
  PRNG prng_elm;
  PRNG prng_mix; /* for mixing function */
  PearsonHash byte_hasher{};
 
 
  using C_t = typename Problem::C_t;
  // using A_t = typename Problem::A_t;
  
  /* Sanity Test:  */
  is_serialize_inverse_of_unserialize<Problem>(Pb, prng_elm);

  
  
  
  /*============================= DICT INIT ==================================*/
  Dict<u64, C_t, Problem> dict{nbytes_memory};
  printf("Initialized a dict with %llu slots = 2^%0.2f slots\n",
	 dict.n_slots, std::log2(dict.n_slots));
  

  /*=============== data extracted from dictionarry ==========================*/
  u64 out0_digest = 0;
  
  
  /*=========================== Collisions counters ==========================*/
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  size_t chain_length1 = 0;

  Counters ctr{}; /* various non-essential counters */
  
  bool found_a_collision = false;
  
  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  bool found_dist = false;


  /* variable to generate families of functions f_i: C -> C */
  /* typename Problem::I_t */
  auto i = Pb.mix_default();


  /* we found a pair of inputs that lead to the golden collisoin or golden claw! */
  bool found_golden_pair = false;

  /*--------------------------debug random -----------------------------------*/
  // /* We need to change to restart calculation with a different function */
  // prng_elm.update_seed(); /* new seed for generatign a mixing function */
  // prng_mix.update_seed();
  // byte_hasher.update_table();
    
  // i = Pb.mix_sample(prng_elm); /* Generates new permutation of f */
  // prng_elm.update_seed(); /* new seed to getting a random value in C_t */
  /* Only for claw: switch the choice between f and g */

  /*----------------------------MAIN COMPUTATION------------------------------*/
  /*=================== Generate Distinguished Points ========================*/
  // todo start here 
  while (not found_golden_pair){

    /* These simulations show that if 10w distinguished points are generated
     * for each version of the function, and theta = 2.25sqrt(w/n) then ...
     */
    /* update F and G by changing `send_C_to_A` and `send_C_to_B` */
    // todo use mix_sample here

    for (size_t n_dist_points = 0;
	 n_dist_points < 10*(dict.n_slots);
	 ++n_dist_points)
      {
      
      found_a_collision = false; /* initially we have not seen anything yet */
      /* fill the input with a fresh random value. */
      Pb.C.randomize(inp_St, prng_elm); /* todo rng should be reviewed */
      
      Pb.C.copy(inp_St, *inp0_pt);
      chain_length0 = 0;

      found_dist = generate_dist_point(Pb,
				       byte_hasher,
				       i, /* permutation number of f's input */
				       chain_length0,
				       difficulty,
				       inp0_pt,
				       out0_pt,
				       inp_mixed,
				       inp0A,
				       args...);/* for claw args := inp0B, inp1B */

      /**************************************************************************/
      #ifdef CLAW_DEBUG
      // GLOBAL VARIABLES FOR DEBUGGING THE NUMBER OF COLLISIONS BEFORE THE GOLDEN
      // 0xdeadbeef tag for debugging
      if (found_golden_A_and_use_f || found_golden_B_and_use_g){
	found_a_collision = true;
	// found_golden_pair = true;
      } else if (Pb.C.is_equal(*out0_pt, Pb.golden_out)){
	/* bad collision, skip it! */
	continue;
      }
      #endif


      #ifdef COLLISION_DEBUG
      // typename Problem::A_t inpA_debug{};
      // Pb.mix(i, *inp0_pt, inp_mixed);
      // Pb.send_C_to_A(inp_mixed, inpA_debug);
      // /* skip adding any to dict except those of the golden inputs */
      // if (not (Pb.is_equal_A(inpA_debug, Pb.golden_inp0)
      // 	       or Pb.is_equal_A(inpA_debug, Pb.golden_inp1)))
      // 	continue;
      // else
      // 	std::cout << "golden value to be inserted!\n"
      // 		  << "inp_unmix = " << *inp0_pt << "\n"
      // 		  << "inpA      = " << inpA_debug << "\n";
      #endif 

      /**************************************************************************/


      out0_digest = Pb.C.hash(*out0_pt);

      
      if (not found_dist)
	continue; /* skip all calculation below and try again  */
      
      ++n_dist_points;
      ctr.increment_n_distinguished_points();
      ctr.n_points += chain_length0;
      
      found_a_collision = dict.pop_insert(out0_digest, /* key */
					   inp_St, /* value  */
					   chain_length0,
					   *inp1_pt,
					   chain_length1,
					   Pb);

      if (found_a_collision) {
	ctr.increment_collisions();
	/* respect the rule that inp0 doesn't have pointers dancing around it */
	Pb.C.copy(inp_St, *inp0_pt); /* (*tmp0_ptO) holds the input value  */

	/*
	 * In case of claw, two extra additional arguments:
	 * (passed as references):
	 * 1) inpA 
	 * 2) inpB
	 */
	found_golden_pair = treat_collision(Pb,
					    byte_hasher,
					    i,
					    inp0_pt,
					    out0_pt,
					    chain_length0,
					    inp1_pt, /* todo fix this */
					    out1_pt,
					    chain_length1,
					    inp_mixed,
					    inp0A, 
					    inp1A,
					    args...); /* for claw args := inp0B, inp1B */
        #ifdef CLAW_DEBUG
	// 0xdeadbeef tag for debugging
	if (found_golden_A_and_use_f && found_golden_B_and_use_g){
	  //found_a_collision = true;
	  found_golden_pair = true;
	}
	#endif


	#ifdef COLLISION_DEBUG
	// if (found_1st_golden_inp and  found_2nd_golden_inp){
	//   std::cout << "Found golden triple in one iteration of the dictioanry\n";
	// }
	  
        #endif 

	
	if (not found_golden_pair)
	    continue; /* nothing to do, test the next one!  */
	else { /* Found the golden pair */
	  print_collision_information(Pb,
				      *out0_pt,
				      *out1_pt,
				      inp0A,
				      inp1A,
				      args... /* inp0B, inp1B*/
				      );


	  double log2_nwords = std::log2(dict.nelements);
          /* todo think about a sensible way to pass |A|, |C|,  */
	  if constexpr (sizeof...(args) == 0) 
	    ctr.save_summary_stats("collision",
				   Pb.nbits_A,/* = |A| */
				   0,/* save summary will ignore this if it 0 */
				   Pb.nbits_C,
				   log2_nwords,
				   difficulty);

	  if constexpr (sizeof...(args) == 2) 
	    ctr.save_summary_stats("claw",
				   Pb.nbits_A,/* = |A| */
				   Pb.nbits_B,/* = |B| */
				   Pb.nbits_C,
				   log2_nwords,
				   difficulty);

	  return; /* nothing more to do */
	}
      }
    }
    /* We need to change to restart calculation with a different function */
    prng_elm.update_seed(); /* new seed for generatign a mixing function */
    prng_mix.update_seed();
    byte_hasher.update_table();
    
    i = Pb.mix_sample(prng_elm); /* Generates new permutation of f */
    prng_elm.update_seed(); /* new seed to getting a random value in C_t */
    /* Only for claw: switch the choice between f and g */

    
    dict.flush();
    ctr.increment_n_updates();

    // 0xdeadbeef tag for debugging
    // reset global variables with each dictionary flush
    #ifdef CLAW_DEBUG
    found_golden_A_and_use_f = false;
    found_golden_B_and_use_g = false;
    #endif

    #ifdef COLLISION_DEBUG
    found_1st_golden_inp = false;
    found_2nd_golden_inp = false;
    #endif 

  }
}

}
#endif
