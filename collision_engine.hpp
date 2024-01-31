#ifndef MITM_COLLISION_ENGINE
#define MITM_COLLISION_ENGINE
#include "AbstractDomain.hpp"
#include "AbstractCollisionProblem.hpp"
#include "base_engine.hpp"

namespace mitm {

  /***************************************************************/
  /* when generating a distinguished point we have:              */
  /*  1)   inp0           =f/g=> out0                            */
  /*  2)  (inp1 := out0)  =f/g=> out1                            */
  /*  3)  (inp2 := out1)  =f/g=> out2                            */
  /*            ...                                              */
  /* m+1) (inp_m := out_m) =f/g=> out_m                          */
  /* A distinguished point found at step `m+1`                   */
  /* "We would like to preserve inp0 at the end of calculation." */
  /***************************************************************/

/*
 * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
 * by out.
 */
template <typename Problem>
void iterate_once(Problem &Pb,
		  typename Problem::C_t& inp,
                  typename Problem::C_t& out,
		  typename Problem::A_t& inpA)
{
        Pb.send_C_to_A(inp, inpA);
        Pb.f(inpA, out);
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
		     typename Problem::A_t& inp0_A,
		     typename Problem::A_t& inp1_A)
{
  /* walk inp0 and inp1 just before `x` */
  /* i.e. iterate_once(inp0) = iterate_once(inp1) */
  /* return false when walking the two inputs don't collide */
  bool found_collision = walk<Problem>(Pb,
				       inp0_chain_len,
				       inp0_pt,
				       out0_pt,
				       inp1_chain_len,
                                       inp1_pt,
				       out1_pt,
				       inp0_A);
  
  /* The two inputs don't lead to the same output */
  if (not found_collision) 
    return false;

  /* useless collision  */
  if ( Pb.C.is_equal(*inp0_pt, *inp1_pt) )
    return false; 

  /* send inputs from C -> A */
  Pb.send_C_to_A(*inp0_pt, inp0_A);
  Pb.send_C_to_A(*inp1_pt, inp1_A);
  
  PAIR_T p{inp0_A, inp1_A};
  collisions_container.push_back(std::move(p)); 
  return true; 
}


/*
 * Given a problem that follows AbstractDomain, and AbstractCollisionProblem.
 * Try to find two inputs x and y s.t. f(x) = f(y), where:
 * f: A -> C
 */
template <typename Problem>
void collisoin_search(Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using C_t = typename Problem::C_t;

  using PAIR_T = std::pair<A_t, A_t>;

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

 
  /* Use these variables to print the full collision */
  A_t inp0A{};
  A_t inp1A{};

  /****************************** EXPERIMENTAL ********************************/

  /* Think about these two varaibles, todo */
  u64  out0_digest = 0; /* hashed value of the output0 */
  /* Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  u8 inp_A_serial[Pb.A.length];
  u8 out0_bytes[Pb.C.length];
  /****************************************************************************/
  int difficulty = 4;
  
  
  /* note the search_engine has different arguments than claw_search */
  search_generic(Pb,
		 collisions_container, /* save found collisions here */
		 1LL<<20, /* #needed_collisions, todo don't hard code it */
	       difficulty,
	       inp0_st, /* starting point in the chain, not a pointer! */
	       inp0_pt,/* pointer to the inp0 s.t. f(inp0) = out0 or using g*/
	       inp1_pt,/* pointer to the 2nd input, s.t. f or g (inp1) = out1 */
	       out0_pt,
	       out1_pt,
	       inp0A,
	       inp1A);

  
}



} /* end of namespace */




#endif
