#ifndef MITM_CLAW_ENGINE
#define MITM_CLAW_ENGINE
#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "engine.hpp"
#include "include/prng.hpp"
#include <boost/stacktrace/stacktrace_fwd.hpp>
#include <stack>
#include <vector>

#include <boost/stacktrace.hpp>
namespace mitm {
/*****************************************************************************/
#ifdef CLAW_DEBUG
template <typename Problem>
void debug_golden_input_A(Problem& Pb,
			  typename Problem::A_t& inpA)
{
  if (Pb.is_equal_A(inpA, Pb.golden_inpA))
    std::cout << "\nwe hit the golden input of A!\n"
		// << boost::stacktrace::stacktrace()
		<< "============================================\n";
}


template <typename Problem>
void debug_golden_input_B(Problem& Pb,
			  typename Problem::B_t& inpB)
{
  if (Pb.is_equal_B(inpB, Pb.golden_inpB))
    std::cout << "\nwe hit the golden input of B!\n"
		// << boost::stacktrace::stacktrace()
		<< "============================================\n";
}


template <typename Problem>
void debug_golden_output(Problem& Pb,
			 typename Problem::C_t& out)
{
  if (Pb.C.is_equal(out, Pb.golden_out))
      std::cout << "\nwe hit the golden output!\n"
		// << boost::stacktrace::stacktrace()
		<< "============================================\n";
}
#endif
/*****************************************************************************/
  
 
template <typename Problem>
void print_collision_information(Problem& Pb,
				 typename Problem::C_t& out0,
				 typename Problem::C_t& out1,
				 typename Problem::A_t& inpA,
				 typename Problem::A_t& inpA_dummy,
				 typename Problem::B_t& inpB_dummy,
				 typename Problem::B_t& inpB)
{
  bool real_collision = Pb.C.is_equal(out0, out1);
  std::cout << "\n++++++++++++++++++++++++++++++++++++++++\n"
	    << "Found golden Pair !\n"
	    << "inpA = " << inpA << "\n"
	    << "inpB = " << inpB << "\n"
	    << "----------------------------------------\n"
	    << "out0 = " << out0 << "\n"
	    << "out1 = " << out1 << "\n"
	    << "out0 == out1? " << real_collision  << "\n"
	    << "++++++++++++++++++++++++++++++++++++++++\n";

}





/* args... = inp0B, inp1B*/
/*
 * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
 * by out_pt.
 */
template <typename Problem>
void iterate_once(Problem &Pb,
		  PearsonHash& byte_hasher,
		  typename Problem::I_t& i, /*permutation number of f's input*/
		  typename Problem::C_t& inp,
                  typename Problem::C_t& out,
		  typename Problem::C_t& inp_mixed,
		  typename Problem::A_t& inpA, /* scratch buffer */
		  typename Problem::B_t& inpB /* scratch buffer */,
		  typename Problem::B_t& inpB_dummy )
{
  /* inpB_dummy only to distinguish this function from iterate_once
   * in collision_engine.hpp. Also, args... will be always passed as two inputs at once
   */

  
  Pb.mix(i, inp, inp_mixed);
  u64 digest = Pb.C.hash(inp_mixed);
  u8 f_or_g = byte_hasher(digest);
  
  if (f_or_g == 1){
    Pb.send_C_to_A(inp_mixed, inpA);
    #ifdef CLAW_DEBUG
    // debug_golden_input_A(Pb, inpA);

    // 0xdeadbeef tag for debugging
    if (Pb.is_equal_A(inpA, Pb.golden_inpA))
      found_golden_A_and_use_f = true;
    //end
    #endif
    
    Pb.f(inpA, out);
  }
  else { /* f_or_g == 0 */
    Pb.send_C_to_B(inp_mixed, inpB);

    #ifdef CLAW_DEBUG
    // debug_golden_input_B(Pb, inpB);

    // 0xdeadbeef tag for debugging 
    if (Pb.is_equal_B(inpB, Pb.golden_inpB))
      found_golden_B_and_use_g = true;
    //end
    #endif
    
    Pb.g(inpB, out);
  }
  #ifdef CLAW_DEBUG
  // debug_golden_output(Pb, out);
  #endif
}




// --------------------------------------------------------------------------------

/*
 * Inputs normall live in A or B, however, the collision only uses C.
 * This function sends the two input to A AND B, if it's not possible, e.g. we can only return 
 */
template <typename Problem>
bool pullback_to_A_B(Problem& Pb,
		     PearsonHash& byte_hasher,
		     typename Problem::I_t& i, /*permutation number of f's input*/
		     typename Problem::C_t& inp0_C,
		     typename Problem::C_t& inp1_C,
		     typename Problem::C_t& inp_mixed,
		     typename Problem::A_t& inp_A,
		     typename Problem::B_t& inp_B)
{
  /* when f is same as g, then there is no point in distinguishingg between */
  /* the two functions. */
  
  /* otherwise, the collision has to be between f and g */
  
  int f_or_g = byte_hasher( Pb.C.hash(inp0_C) );
  
  if (f_or_g == 1){ 
    if (byte_hasher( Pb.C.hash(inp1_C) )  == 0){
      /* inp0 -> A, inp1 -> B */
      Pb.mix(i, inp0_C, inp_mixed);
      Pb.send_C_to_A(inp_mixed, inp_A);

      Pb.mix(i, inp1_C, inp_mixed);
      Pb.send_C_to_B(inp_mixed, inp_B);
    }
  } else { 
    if ( byte_hasher( Pb.C.hash(inp1_C) ) == 1){ /* a sensible case*/
      /* inp0 -> B, inp1 -> A */
      Pb.mix(i, inp1_C, inp_mixed);
      Pb.send_C_to_A(inp_mixed, inp_A);

      Pb.mix(i, inp0_C, inp_mixed);
      Pb.send_C_to_B(inp_mixed, inp_B);
      return true;
    }
  }
  return false;
}


  
/*
 * Given a potential collision, walk the two inputs until reaching a common
 * point. If it's a claw problem, check that the two walked inputs reach the
 * same point using two different functions, f and g. Also, check it is not a
 * robinhood, and use the provided golden collision criterion.
 */ 
template <typename Problem>  
bool treat_claw(Problem& Pb,
		     PearsonHash& byte_hasher,
		     typename Problem::I_t& i,
		     typename Problem::C_t*& inp0_pt,
		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::C_t*& inp1_pt,
		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
		     typename Problem::C_t& inp_mixed,
		     typename Problem::A_t& inp0_A,
		     typename Problem::A_t& inp1_A, /* <- dummy arg */  
		     typename Problem::B_t& inp0_B, /* <- dummy arg */ 
		     typename Problem::B_t& inp1_B)
{

  /*
   * Dummy inputs distinguishes make claw's treat_collision always distinguish 
   * from collisions' treat_collision, even when the types A_t = B_t
   */
  /* walk inp0 and inp1 just before `x` */
  /* i.e. iterate_once(inp0) = iterate_once(inp1) */
  /* return false when walking the two inputs don't collide */
  bool found_collision = walk<Problem>(Pb,
				       byte_hasher,
				       i,
				       inp0_chain_len,
				       inp0_pt,
				       out0_pt,
				       inp1_chain_len,
                                       inp1_pt,
				       out1_pt,
				       inp_mixed,
				       inp0_A,
				       inp0_B, // this arg and one below allows
				       inp1_B); // us to call walk for claw.

  /* The two inputs don't lead to the same output */
  if (not found_collision) 
    return false;

  /* useless collision  */
  if ( Pb.C.is_equal(*inp0_pt, *inp1_pt) )
    return false;

  /*  f: A -> C, g: B -> C. Can we pullback the inputs *inp0 and *inp1_pt in C
   *  to inp_A in A and inp_B (the order doesn't matter)? If yes, write the
   * results to inp0_A and inp0_B. Otherwise, return false
   */
  bool is_potential_collision = pullback_to_A_B(Pb,
						byte_hasher,
						i,
						*inp0_pt, // a given input in C 
						*inp1_pt, // a given inpunt in C
						inp_mixed,
						inp0_A, // resulted input in A
						inp1_B);// resulted input in B
  if (not is_potential_collision)
    return false; /* don't add this pair */
  
  /* */
  return Pb.is_good_pair(*out0_pt, inp0_A, inp1_B);
}


}

#endif
