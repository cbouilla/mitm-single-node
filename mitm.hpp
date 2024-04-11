#ifndef MITM
#define MITM

#include "AbstractDomain.hpp"
#include "AbstractClawProblem.hpp"
#include "AbstractCollisionProblem.hpp"
#include "engine.hpp"
#include "claw_engine.hpp"
#include "collision_engine.hpp"
#include "naive_engine.hpp"


/* Interface for the search functions:
 * - collision_search
 * - claw_search
 * Or the naive search which is guaranteed to find all collisions/claws
 * 
 * 
 */
namespace mitm {


/*
 * Given a problem that follows AbstractDomain, and AbstractCollisionProblem.
 * Try to find two inputs x and y s.t. f(x) = f(y), where:
 * f: A -> C
 */
template <typename Problem>
void collision_search(Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using C_t = typename Problem::C_t;


  /* ============================= BUFFERS ================================== */
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
  A_t inp1A{};

  /****************************** EXPERIMENTAL ********************************/

  /* Think about these two varaibles, todo */
  // u64  out0_digest = 0; /* hashed value of the output0 */
  /* Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  // u8 inp_A_serial[Pb.A.length];
  // u8 out0_bytes[Pb.C.length];
  /****************************************************************************/
  int difficulty = 0;
  
  
  /* note the search_engine has different arguments than claw_search */
  search_generic(Pb,
		 difficulty,
		 inp0_st, /* starting point in the chain, not a pointer! */
		 inp0_pt,/* pointer to the inp0 s.t. f(inp0) = out0 or using g*/
		 inp1_pt,/* pointer to the 2nd input, s.t. f or g (inp1) = out1 */
		 out0_pt,
		 out1_pt,
		 inp_mixed,
		 inp0A,
		 inp1A);

  
}




/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 */
template <typename Problem>
void claw_search(Problem& Pb, int difficulty = 0)
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;


  /* ============================= BUFFERS ================================== */

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

  
  search_generic(Pb,
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


/*
 * Given a problem that follows AbstractDomain, and AbstractCollisionProblem.
 * Try to find two inputs x and y s.t. f(x) = f(y), where:
 * f: A -> C using a naive algortihm
 * Note: It's guaranteed in this method to find all collisions.
 * This function asks for an extra input "ith_element"
 * s.t. ith_element(x, i), x <- ith_input (according to some order, doesn't matter which )
 *      ith_element(x, i) =/= ith_element(x, j) if i =/= j.
 */
template <typename Problem>
auto naive_collisoin_search(Problem& Pb,
			    size_t n_elements
			    ) -> std::vector<std::pair<typename Problem::C_t,
						       typename Problem::A_t>>
{

  using A_t = typename Problem::A_t;
  using C_t = typename Problem::C_t;
  std::function<void(const A_t&, C_t&)> f = [&Pb](const A_t& x, C_t& y) { Pb.f(x, y); };
  std::function<void(A_t const&, size_t)> ith_element = [&Pb](const A_t& x, size_t i) { Pb.ith_element(x, i); };
  
  /* create all possible pairs (f(x), x) sorted by f(x) as
   * std::vector<std::pair<>>
   */

  
  auto inps_outs = domain_images<A_t, C_t>(f,
					   ith_element,
					   0, /* start, this argument for future use with MPI */
					   n_elements);

  /* return all pairs (f(x), x) s.t. there is an x' =/= x and f(x') == f(x), do this for all x */
  auto collisions = extract_collisions(inps_outs, inps_outs);

  /* save collisions in disk */
  // todo
  return collisions;
}


/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 * Note: It's guaranteed in this method to find all claws.
 */
template <typename Problem>
auto naive_claw_search(Problem &Pb,
		       size_t n_elements_A,
		       size_t n_elements_B
		       ) -> std::vector<typename Problem::C_t>
{

  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;

  /* since there is this argument to f, and g */
  std::function<void(const A_t&, C_t&)> f = [&Pb](const A_t& x, C_t& y) { Pb.f(x, y); };
  std::function<void(const B_t&, C_t&)> g = [&Pb](const B_t& x, C_t& y) { Pb.g(x, y); };
  std::function<void(A_t &, size_t)> ith_element = [&Pb](const A_t& x, size_t i) { Pb.ith_element(x, i); };

  
  /* get all f(x) where x in A, stored in std::vector<C_t>  */
  auto f_images = images(f,
			 ith_element,
			 0, /* start with the 1st element in A */
			 n_elements_A); /* end with the last element in A */

  /* get all f(x) where x in A, stored in std::vector<C_t>  */
  auto g_images = images(g,
			 ith_element,
			 0, /* start with the 1st element in A */
			 n_elements_B); /* end with the last element in A */
  
  /* todo iterating only over A inputs, may not return all claws!
   *  we should iterate over the largest domain
   */
  auto collisions = extract_collisions(f_images,
				       g_images);

  /* todo save the results in disk! */


  return collisions;
  
}

}

#endif
