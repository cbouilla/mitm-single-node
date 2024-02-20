#ifndef MITM_CLAW_ENGINE
#define MITM_CLAW_ENGINE
#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "base_engine.hpp"
#include <vector>


namespace mitm {

/* args... = inp0B, inp1B*/
/*
 * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
 * by out_pt.
 */
template <typename Problem>
void iterate_once(Problem &Pb,
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
  int f_or_g = Pb.C.extract_1_bit(inp_mixed);

  if (f_or_g == 1){
    Pb.send_C_to_A(inp_mixed, inpA);
    Pb.f(inpA, out);
  }
  else { /* f_or_g == 0 */
    Pb.send_C_to_B(inp_mixed, inpB);
    Pb.g(inpB, out);
  }
}





// --------------------------------------------------------------------------------

/*
 * Inputs normall live in A or B, however, the collision only uses C.
 * This function sends the two input to A AND B, if it's not possible, e.g. we can only return 
 */
template <typename Problem>
bool send_2_A_and_B(Problem& Pb,
		    typename Problem::C_t& inp0_C,
		    typename Problem::C_t& inp1_C,
		    typename Problem::A_t& inp_A,
		    typename Problem::B_t& inp_B)
{
  /* when f is same as g, then there is no point in distinguishingg between */
  /* the two functions. */
  
  /* otherwise, the collision has to be between f and g */
  int f_or_g = Pb.C.extract_1_bit(inp0_C);
  if (f_or_g == 1){ 
    if (Pb.C.extract_1_bit(inp1_C) == 0){ /* a sensible case*/
      Pb.send_C_to_A(inp0_C, inp_A);
      Pb.send_C_to_B(inp1_C, inp_B);
    }
  } else { /* Pb.C.extract_1_bit(inp0_C) == 0*/
    if (Pb.C.extract_1_bit(inp1_C) == 1){ /* a sensible case*/
      // good case
      Pb.send_C_to_A(inp1_C, inp_A);
      Pb.send_C_to_B(inp0_C, inp_B);
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
 */ // todo: the pair type looks ugly, rethink the solution.
template <typename Problem, typename PAIR_T> /* PAIR_T = std::pair<A_t, B_t> */
bool treat_collision(Problem& Pb,
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
				       i,
				       inp0_chain_len,
				       inp0_pt,
				       out0_pt,
				       inp1_chain_len,
                                       inp1_pt,
				       out1_pt,
				       inp_mixed,
				       inp0_A,
				       inp0_B,
				       inp1_B);

  /* The two inputs don't lead to the same output */
  if (not found_collision) 
    return false;

  /* If found a robinhood, don't  */
  if ( Pb.C.is_equal(*inp0_pt, *inp1_pt) )
    return false;

  /* when f=/= g the one of the input should an A input while the other is B's */
  /* If this is not satisfied return */
  bool is_potential_collision = send_2_A_and_B(Pb,
					       *inp0_pt,
					       *inp1_pt,
					       inp0_A,
					       inp1_B);
  if (not is_potential_collision)
    return false; /* don't add this pair */
  
  
  return Pb.is_good_pair(*out0_pt, inp0_A, inp0_B);
}



/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 */
template <typename Problem>
void claw_search(Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;

  using PAIR_T = std::pair<A_t, B_t>;

  /* ============================= BUFFERS ================================== */
  std::vector<PAIR_T> collisions_container{};


  /* Input/Output containers */
  /* 1st set of buffers: Related to input0 as a starting point */
  /* either tmp0 or  output0 */
  C_t inp0_or_out0_buffer0{};
  C_t inp0_or_out0_buffer1{};


  /* 2nd set of buffers: Related to input1 as a starting point */
  /* When we potentially find a collision, we need 2 buffers for (inp1, out1) */
  /* We will use the initial value of inp1 once, thus we don't need to gaurd  */
  /* in an another variable untouched */
  C_t inp1_or_out1_buffer0{};
  C_t inp1_or_out1_buffer1{};


  /* -------------------------- Passed variables ---------------------------- */
  C_t inp0_st{}; /* starting input on the chain */

  /* Always points to the region that contains the output,  */
  C_t* inp0_pt = &inp0_or_out0_buffer0; /* even if it's not the same address*/
  C_t* inp1_pt = &inp1_or_out1_buffer0;

  /* Always points to the region that contains the output */
  C_t* out0_pt = &inp0_or_out0_buffer1;/* even if it's not the same address*/
  C_t* out1_pt = &inp1_or_out1_buffer1;

 
  C_t inp_mixed{};
  
  
  /* Use these variables to print the full collision */
  A_t inp0A{};
  B_t inp0B{};
  A_t inp1A{};
  B_t inp1B{};

  /****************************** EXPERIMENTAL ********************************/
  /* Think about these two varaibles, todo */
  // u64  out0_digest = 0; /* hashed value of the output0 */
  /* Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  // u8 inp_A_serial[Pb.A.length];
  // u8 inp_B_serial[Pb.B.length];
  // u8 out0_bytes[Pb.C.length];
  /****************************************************************************/
  int difficulty = 4;
  
  search_generic(Pb,
		 collisions_container, /* save found collisions here */
		 1LL<<20, /* #needed_collisions, todo don't hard code it */
		 difficulty,
		 inp0_st, /* starting point in the chain, not a pointer! */
		 inp0_pt,/* pointer to the inp0 s.t. f(inp0) = out0 or using g*/
		 inp1_pt,/* pointer to the 2nd input, s.t. f or g (inp1) = out1 */
		 out0_pt,
		 out1_pt,
		 inp_mixed,
		 inp0A,
		 inp1A,
		 inp0B, /* Last two inputs are args... in search generic */
		 inp1B);

  
}

}

#endif
