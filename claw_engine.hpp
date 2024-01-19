#ifndef MITM_CLAW_ENGINE
#define MITM_CLAW_ENGINE

namespace mitm {
// for now keep this code

/*
 * Inputs normall live in A or B, however, the collision only uses C.
 * This function sends the two input to A AND B, if it's not possible, e.g. we can only return 
 */
template <typename Problem>
bool send_2_A_and_B(typename Problem::C_t& inp0_C,
		    typename Problem::C_t& inp1_C,
		    typename Problem::A_t& inp_A,
		    typename Problem::B_t& inp_B,
		    Problem& Pb)
{
  /* when f is same as g, then there is no point in distinguishingg between */
  /* the two functions. */
  if constexpr(Pb.f_eq_g == 1){
    Pb.send_C_to_A(inp0_C, inp_A);
    Pb.send_C_to_B(inp1_C, inp_B);
    return true;
  }
  
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
 * return false if the two inputs leads to a robinhood or collision on the same
 * function (when f =/= g). todo: walk should also return false when it could
 * not get a collisison
 */
template <typename Problem >
bool treat_collision(typename Problem::C_t*& inp0_pt,
		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::C_t*& inp1_pt,
		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
                     std::vector< std::pair<typename Problem::A_t,
		                            typename Problem::B_t> >& container,
		     Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;

  /* walk inp0 and inp1 just before `x` */
  /* i.e. iterate_once(inp0) = iterate_once(inp1) */
  /* return false when walking the two inputs don't collide */
  bool found_collision = walk<Problem>(inp0_pt,
				       out0_pt,
				       inp0_chain_len,
				       inp1_pt,
				       out1_pt,
				       inp1_chain_len,
				       Pb);

  /* The two inputs don't lead to the same output */
  if (not found_collision) 
    return false;

  /* If found a robinhood, don't  */
  if ( Pb.C.is_equal(*inp0_pt, *inp1_pt) )
    return false; 

  /* send inp{0,1} to inp_A and inp_B. If it is not possible, return false.*/
  A_t inp_A{};
  B_t inp_B{};


  /* when f=/= g the one of the input should an A input while the other is B's */
  /* If this is not satisfied return */
  bool is_potential_collision = send_2_A_and_B(*inp0_pt, *inp1_pt,
					       inp_A,     inp_B,
					       Pb);
  if (not is_potential_collision)
    return false; /* don't add this pair */
  
  std::pair<A_t, B_t> p{inp_A, inp_B};
  container.push_back(std::move(p)); 
  return true;
}

  
}

#endif
