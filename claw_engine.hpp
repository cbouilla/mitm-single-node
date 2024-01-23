#ifndef MITM_CLAW_ENGINE
#define MITM_CLAW_ENGINE
#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "engine.hpp"


namespace mitm {

template <typename Problem>
void iterate_once(Problem &Pb,
		  typename Problem::C_t& inp,
                  typename Problem::C_t& out,
		  typename Problem::A_t& inpA, /* scratch buffer */
		  typename Problem::B_t& inpB /* scratch buffer */ )
{
  /*
   * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
   * by out_pt.
   */
  int f_or_g = Pb.C.extract_1_bit(inp);
  if (f_or_g == 1){
    Pb.send_C_to_A(inp, inpA);
    Pb.f(inpA, out);
  }
  else { /* f_or_g == 0 */
    Pb.send_C_to_B(inp, inpB);
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
		     std::vector<PAIR_T> collisions_container,
		     typename Problem::C_t*& inp0_pt,
		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::C_t*& inp1_pt,
		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
		     typename Problem::A_t*& inp0_A,
		     typename Problem::A_t*& inp1_A, /* dummy input */
		     typename Problem::B_t*& inp0_B, /* dummy input */
		     typename Problem::B_t*& inp1_B,
		     int dummy_arg /* To distinguish it from claw treat... */)
{

  /*
   * Dummy inputs distinguishes make claw's treat_collision always distinguish 
   * from collisions' treat_collision, even when the types A_t = B_t
   */
  /* walk inp0 and inp1 just before `x` */
  /* i.e. iterate_once(inp0) = iterate_once(inp1) */
  /* return false when walking the two inputs don't collide */
  bool found_collision = walk<Problem>(Pb,
				       inp0_chain_len,
				       inp0_pt,
				       out0_pt,
				       inp1_chain_len,
                                       inp1_pt,
				       out1_pt);


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
					       *inp0_A,
					       *inp1_B);
  if (not is_potential_collision)
    return false; /* don't add this pair */
  
  PAIR_T p{inp0_A, inp1_B};
  collisions_container.push_back(std::move(p)); 
  return true;
}



/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 */
template <typename Problem>
void claw_search(Problem& Pb)
{
  
}

  
}

#endif
