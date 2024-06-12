#ifndef MITM_MPI_ENGINE
#define MITM_MPI_ENGINE
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
#include "sender.hpp"
#include "receiver.hpp"


namespace mitm {


/* Parallel search finds collision/claw of
 * iterate_once: C -> C, these collisions can be converted into a solution
 * to the problem we are treating, claw or a collision.
 *  
 */
template<typename Problem,  typename... Types>
void parallel_search(Problem& Pb,
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
  //                               digest        chain_length 
  size_t triple_size = Pb.C.size + sizeof(u64) + sizeof(u32);
  int nreceivers = get_nodes_count(); // for now make them = nservers 
  int difficulty = 0; // todo Automatically find the best difficulty.

  MITM_MPI_data my_info = MITM_MPI_Init(nreceivers,
					triple_size, //  |(inp, hash, chain_length)|
					0);

  if (my_info.my_role == SENDER)
    sender(Pb,
	   my_info,
	   difficulty,
	   inp_St,
	   inp0_pt,
	   inp1_pt,
	   out0_pt,
	   out1_pt,
	   inp_mixed,
	   inp0A,
	   inp1A,
	   args...);
  
  if (my_info.my_role == RECEIVER)
    reciever(Pb,
	     my_info,
	     difficulty,
	     inp_St,
	     inp0_pt,
	     inp1_pt,
	     out0_pt,
	     out1_pt,
	     inp_mixed,
	     inp0A,
	     inp1A,
	     args...);
}
}
#endif 
