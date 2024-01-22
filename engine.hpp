#ifndef MITM_ENGINE
#define MITM_ENGINE



#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "AbstractCollisionProblem.hpp"
#include "include/prng.hpp"
#include "include/timing.hpp"
#include "include/memory.hpp"
#include "dict.hpp"

#include <exception>
#include <iostream>
#include <vector>

/// When a function, func, behavior's depends if the problem is a claw or a
///  collision, we use veriadic template to have two different implementation:
/// template<_, typename... Types>  /* A veriadice template */
/// func(_, Types... args)
/// Then in `claw_engine.hpp` would have a template specialization of func
/// with more arguments than the specialization of func found in
/// `collisions_engine.hpp`
/// List of functions that are different for claw and collision problem:
/// 1) `iterate_once` claw_engine.hpp introduces 2 more arguments.
/// 

// breaking code notes:
// 1- the order of arguments have changed. Pb is the first now.
// 2- inp0 should be copied outside the function.

namespace mitm {


/*
 * Defines how to walk in a sequence. Use functions overloading to differentiate
 * between a claw's iteration and a collisions's iteration. Specfically:
 * When it's a claw problem, then it would have 2 extra arguments, namely:
 * 1) inp_A 2) inp_B. Besides Pb, inp, and out.
 * When it's a collision, args would be empty. i.e. it only uses:
 * Pb, inp, and out.
 */
template <typename Problem, typename... Types>
void iterate_once(Problem &Pb,
		  typename Problem::C_t &inp,
                  typename Problem::C_t &out,
		  Types... args)
{
  /* C++ always prefers more specialized templates. 
   * For the claw code, see `iterate_once` in `claw_engine.hpp`
   * For the collsions code, see `iterate_once` in `collision_engine.hpp`
   */
  std::cerr << "Internal Error: should not use general implementation of `iterate_once`!\n";
  std::terminate(); /* Never use this implementation! */
}



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
			 u64& chain_length, /* write chain lenght here */
			 const i64 difficulty, /* #bits are zero at the beginning */
			 typename Problem::C_t*& inp_pt,
			 typename Problem::C_t*& out_pt,
			 Types... args)
{
  /* todo break note copy inp0 before passing it here! */
  // code below was uncommented.
  // /* copy the input to tmp, then never touch the inp again! */
  // Pb.C.copy(inp0, *inp_pt);

  const u64 mask = (1LL<<difficulty) - 1;
  u64 digest = 0;
  bool found_distinguished = false;
  
  /* The probability, p, of NOT finding a distinguished point after the loop is
   * Let: theta := 2^-d
   * difficulty, N = k*2^difficulty then,
   * p = (1 - theta)^N =>  let ln(p) <= -k
   */
  constexpr u64 k = 40;
  for (u64 i = 0; i < k*(1LL<<difficulty); ++i){
    /* uses claw's iterate_once if args... is not empty, otherwise collisions'*/
    iterate_once(Pb, *inp_pt, *out_pt, args...); 
    ++chain_length;

    /* we may get a dist point here */
    digest = Pb.C.hash(*out_pt);
    found_distinguished = is_distinguished_point(digest, mask);

    /* no need to continue if a distinguished is found  */
    if (found_distinguished) [[unlikely]]{ return true; }

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
      /* todo throw exception */
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
template<typename Problem>
bool walk(Problem& Pb,
	  u64 inp0_chain_len,
	  typename Problem::C_t*& inp0_pt,
	  typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
	  u64 inp1_chain_len,
          typename Problem::C_t*& inp1_pt,
	  typename Problem::C_t*& out1_pt /* inp1 calculation buffer */)
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
    iterate_once(*inp0_pt, *out0_pt, Pb);
    swap_pointers(inp0_pt, out0_pt);
  }
  
  for (; inp0_chain_len < inp1_chain_len; --inp1_chain_len){
    iterate_once(*inp1_pt, *out1_pt, Pb);
    swap_pointers(inp1_pt, out1_pt);
  }

  /*****************************************************************************/
  /* now both inputs have equal amount of steps to reach a distinguished point */
  /* both sequences needs exactly `len` steps to reach distinguished point.    */
  for (size_t i = 0; i < len; ++i){
    /* walk them together and check each time if their output are equal     */
    /* return as soon equality is found. The equality could be a robinhood. */
    iterate_once(*inp0_pt, *out0_pt, Pb);
    iterate_once(*inp1_pt, *out1_pt, Pb);
    
    if(Pb.C.is_equal( *out0_pt, *out1_pt ))
      return true; /* They are equal */

    swap_pointers(inp0_pt, out0_pt);
    swap_pointers(inp1_pt, out1_pt);
  }
  return false;
}


/*
 * In case of claw: two extra additional arguments:
 * 1) inpA_pt
 * 2) inpB_pt
 */
template <typename Problem, typename PAIR_T, typename... Types >
bool treat_collision(Problem& Pb,
		     std::vector<PAIR_T> collisions_container,
		     typename Problem::C_t*& inp0_pt,
		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::C_t*& inp1_pt,
		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
		     Types... args)
{
  std::cerr << "Internal Error: should not use general implementation of if `treate_collision`!\n";
  std::terminate(); /* Never use this implementation! */
}

// Start from here

template<typename Problem, typename PAIR_T, typename... Types>
PAIR_T search_generic(Problem& Pb,
		      typename Problem::C_t& inp_St, // Startign point in chain
		      typename Problem::C_t* inp0_pt,
		      typename Problem::C_t* inp1_pt,
		      typename Problem::C_t* out0_pt,
		      typename Problem::C_t* out1_pt,
		      typename Problem::A_t* inp0A_pt,
		      typename Problem::A_t* inp1A_pt,
		      Types... args) /* two extra arguments if claw problem
				      * 1) inp0B_pt
				      * 2) inp1B_pt
				      */
{
  PRNG rng_urandom;

  using C_t = typename Problem::C_t;
  using A_t = typename Problem::A_t;
  
  /* Sanity Test:  */
  is_serialize_inverse_of_unserialize<Problem>(Pb, rng_urandom);

  /*----------------------------- DICT INIT ----------------------------------*/
  size_t n_bytes = 0.5*get_available_memory();
  std::cout << "Going to use "
	    << std::dec << n_bytes << " bytes = 2^"<< std::log2(n_bytes)
	    << " bytes for dictionary!\n";
  

  Dict<u64, C_t, Problem> dict{n_bytes};
  /* In claw: PAIR_T = <A_t, B_t>, in collisions PAIR_T = <C_t, C_t> */
  std::vector<PAIR_T> collisions_container{};
  u64 out0_digest = 0;
  
  
  std::cout << "Initialized a dict with " << dict.n_slots << " slots = 2^"
	    << std::log2(dict.n_slots) << " slots\n";

  /*--------------------------- Collisions counters --------------------------*/

  int difficulty = 9;
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  size_t chain_length1 = 0;
  
  bool is_collision_found = false;
  size_t n_collisions = 0;
  size_t n_needed_collisions = 1LL<<20;
  
  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  size_t n_robinhoods = 0;
  bool found_dist = false;
  size_t n_distinguished_points = 0;
  constexpr size_t interval = (1LL<<15);
  double collision_timer = wtime();
  

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
      Pb.C.randomize(inp_St, rng_urandom);
      Pb.C.copy(inp_St, *inp0A_pt);
      chain_length0 = 0; 
      
      
      found_dist = generate_dist_point<Problem>(inp0_pt,
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
					   inp_St, /* value  */
					   chain_length0,
					   *inp1_pt, /* save popped element here */
					   chain_length1, /* of popped input */
					   Pb);
      
      if (is_collision_found) [[unlikely]]{
	++n_collisions;

	/* Move this code to print collision information */
        std::cout << "\nA collision is found\n"
		  << "It took " << (wtime() - collision_timer) << " sec\n"
		  << "inp0 (starting point) = " << inp_St << "\n"
		  << "digest0 = 0x" << out0_digest << "\n"
		  << "chain length0 = " << chain_length0 << "\n"
		  << "inp1 (starting point) = " << *inp1_pt << "\n"
	  	  << "chain length1 = " << std::dec << chain_length1 << "\n"
		  << "-------\n";

	collision_timer = wtime();
	/* respect the rule that inp0 doesn't have pointers dancing around it */
	Pb.C.copy(inp_St, *inp0_pt); /* (*tmp0_ptO) holds the input value  */

	/*
	 * In case of claw: two extra additional arguments:
	 * 1) inpA_pt
	 * 2) inpB_pt
	 */
	bool is_potential_coll = treat_collision<Problem>(Pb,
							  collisions_container,
							  inp0_pt,
							  out0_pt,
							  chain_length0,
							  inp1_pt,
							  out1_pt,
							  chain_length1,
							  args...);

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
	}
      }
    }
  }
  /* end of work */
  /* todo fix this return type! */
  return PAIR_T(*inp0_pt, *inp1_pt); // todo wrong values
}


}

#endif
