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
enum MITM_MPI_TAGS {INTERCOM_TAG}; // to be extended ...


struct MITM_MPI_data{
  MPI_Comm global_comm; // Global communicator. Kept just in case, todo to be removed!
  MPI_Comm local_comm;  // Processes that do the same thing, e.g. senders. 
  MPI_Comm inter_comm;  // splitted into two: senders and receivers.
  size_t const nelements_buffer = 1000; // #elements stored in buffer during send/receive
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
  
  /*  Receivers establish their communication channel with senders */
  if (is_receiver) {
    MPI_Intercomm_create(my_info.local_comm, /* my tribe */
			 0, /* Local name of the sheikh of my tribe */
			 MPI_COMM_WORLD, /* Where to find the remote leader. */
			 nreceivers, /* Global Name/rank of senders leader. */
			 INTERCOM_TAG, /* This tag is unique to the creation of intercomm */
			 &my_info.inter_comm); /* <- put down these info here */

    my_info.my_role = RECEIVER;
  }
  
  /* senders establish their communication channel with receivers */
  else {
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
  


  return my_info;
}


/* Return how many bytes of RAM can a receiver allocate,
 * if the available space is shared. */
size_t get_nbytes_per_receivers(MITM_MPI_data& my_info)
{
  return 0;
}

/* One round of computation for sender */
template <typename Problem, typename... Types> void sender_round() {}

/* One round of computation of receiver  */
template <typename Problem, typename... Types> void receiver_round() {}


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


  /* ----------------------------- sender buffer ---------------------------- */
  u8* snd_buff;
  snd_buff = new u8[my_info.nreceivers * my_info.nelements_buffer * Pb.C.size ];

  // ======================================================================== //
  //                                 STEP 1                                   //
  // ------------------------------------------------------------------------ //
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

   while (true){
    sender_round<typename Problem, typename Types>();

    /* Update seeds  */

    /* Clear buffers */
    std::memset(snd_buff, 0, my_info.nreceivers * my_info.nelements_buffer * Pb.C.size);
    /* repeat! */
  }
   delete[] snd_buff;
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
  // ======================================================================== //
  //                          STEP 0:  INTITALIZATION                         //
  // ------------------------------------------------------------------------ //

  
  /* ------------------------------- DICT INIT ------------------------------ */
  size_t nbytes_memory = get_nbytes_per_receivers(my_info);
  
  Dict<u64, C_t, Problem> dict{nbytes_memory};
  printf("Recv %d initialized a dict with %llu slots = 2^%0.2f slots\n",
	 my_info.my_rank_intercomm, dict.n_slots, std::log2(dict.n_slots));

  /* ----------------------------- receive buffer ---------------------------- */
  u8* rcv_buff;
  rcv_buff = new u8[my_info.nelements_buffer * Pb.C.size];

  /* ---------------- data inserted/extracted from dictionarry -------------- */
  u64 out0_digest = 0;
  size_t chain_length1 = 0;
  // C_t* out1_pt,
  
  /*--------------------------- Collisions counters --------------------------*/
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  Counters ctr{}; /* various non-essential counters */

  /* ----------------------------- Query Results ---------------------------- */
  bool found_a_collision = false;
  
  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  bool found_dist = false;

  /* we found a pair of inputs that lead to the golden collisoin or golden claw! */
  bool found_golden_pair = false;



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

  while (true){
    receiver_round<typename Problem, typename Types>();

    /* Update seeds  */
    /* Clear receiv buffers */
    std::memset(rcv_buff, 0, my_info.nelements_buffer * Pb.C.size);
    /* Clear dictionary */
    dict.flush();
    /* repeat! */
  }

  delete[] rcv_buff;
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


	  double log2_nbytes = std::log2(nbytes_memory);
          /* todo think about a sensible way to pass |A|, |C|,  */
	  if constexpr (sizeof...(args) == 0) 
	    ctr.save_summary_stats("collision",
				   Pb.nbits_A,/* = |A| */
				   0,/* save summary will ignore this if it 0 */
				   Pb.nbits_C,
				   log2_nbytes,
				   difficulty);

	  if constexpr (sizeof...(args) == 2) 
	    ctr.save_summary_stats("claw",
				   Pb.nbits_A,/* = |A| */
				   Pb.nbits_B,/* = |B| */
				   Pb.nbits_C,
				   log2_nbytes,
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
