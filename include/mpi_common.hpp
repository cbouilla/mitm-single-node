#ifndef MITM_MPI_COMMON
#define MITM_MPI_COMMON

#include <cstddef>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <mpi.h>

#include "../docs/AbstractDomain.hpp"
#include "../docs/AbstractClawProblem.hpp"
#include "../docs/AbstractCollisionProblem.hpp"

#include "util/prng.hpp"
#include "util/memory.hpp"

#include "counters.hpp"
#include "dict.hpp"
#include "engine.hpp"

/******************************************************************************/
/*                                     MPI                                    */
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


namespace mitm {
enum process_type{ SENDER,  RECEIVER };
// enum MITM_MPI_TAGS {INTERCOM_TAG, ROUND_SND_TAG}; // to be extended ...
#define INTERCOM_TAG 0
#define ROUND_SND_TAG 1
  
/* A wrapper to MPI data associated with a process. */
struct MITM_MPI_data{
  MPI_Comm global_comm; // Global communicator. Kept just in case, todo to be removed!
  MPI_Comm local_comm;  // Processes that do the same thing, e.g. senders. 
  MPI_Comm inter_comm;  // splitted into two: senders and receivers.
  size_t const nelements_buffer = 1000; // #elements stored in buffer during send/receive
  size_t nelements_mem_receiver; // How many elements can be stored across all dictionaries.
  size_t msg_size = 0; // todo initialize this variable. // Size of msg with tag ROUND_SND_TAG + i where i = 0, 1, ...
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
 * It requires MPI-3 or higher.
 * source: https://stackoverflow.com/a/34118174 */
int get_nodes_count()
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

/* Return how many bytes a sender would allocated for sending buffers */
size_t sender_buffer_size(MITM_MPI_data& my_info,
			  size_t element_size /* |(inp, digest, chain_length)| */
			  )
{
  return (my_info.nreceivers + 1) * (my_info.nelements_buffer)
         * element_size + sizeof(u8);
}
  
/* Initializes mpi processes for mitm,  and splits senders and receivers into
 * seperate communicators. Aslo, gives the number of total numbers of senders
 *  and receivers, `nsenders` and `nreceivers` respectively,  and the number
 * senders and receivers that are on the same node, `nsenders_node` & 
 * `nreceivers_node` respectively. This allows us to use nodes with different
 * specs without affecting the code.
 * 
 */
MITM_MPI_data MITM_MPI_Init(int nreceivers,
			    size_t element_size, //  |(inp, hash, chain_length)|
			    int argc=0, char** argv = NULL)
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
		 ); /* sad face for the separation between processes */


  /* Create communicator between the two clans of processes. */
  if (is_receiver) {
    /*  Receivers establish their communication channel with senders */
    MPI_Intercomm_create(my_info.local_comm, /* my tribe */
			 0, /* Local name of the sheikh of my tribe */
			 MPI_COMM_WORLD, /* Where to find the remote leader. */
			 nreceivers, /* Global rank of senders leader. */
			 INTERCOM_TAG, /* This tag is unique to the creation of intercomm */
			 &my_info.inter_comm); /* <- put down these info here */
    my_info.my_role = RECEIVER;
    
  } else {
    /* senders establish their communication channel with receivers */
    MPI_Intercomm_create(my_info.local_comm, /* my tribe */
			 0, /* Local name of the leader of my tribe */
			 MPI_COMM_WORLD, /* Where to find the remote leade. */
			 0, /* Global rank of receivers leader. */
			 INTERCOM_TAG, /* This tag is unique to the creation of intercomm */
			 &my_info.inter_comm); /* <- put down these info here */

    my_info.my_role = SENDER;
  }

  /* Fill the rest of the details  */
  my_info.nsenders = my_info.nprocesses - nreceivers;
  my_info.nreceivers = nreceivers;
  my_info.nnodes = get_nodes_count();
  my_info.nsenders_node = my_info.nsenders / my_info.nnodes;
  my_info.nreceivers_node = my_info.nreceivers / my_info.nnodes;
  
  MPI_Comm_rank(my_info.local_comm, &my_info.my_rank_intercomm);
  

  /* We can say how many elements can be stored in memory */
  // todo formulas are incorrect!
  size_t nbytes_avaiable_node = get_available_memory();
  /* How many bytes a sender needed for the sending buffers */
  /* On my local node, how many bytes senders would consume? */
  size_t senders_needed_bytes_node = my_info.nsenders_node
                                   * sender_buffer_size(my_info, element_size);
  

  //   size_t dict_one_element_size = element_size;
  // memory available to all receivers per node 
  size_t nelements_ram_receivers_node = (nbytes_avaiable_node - senders_needed_bytes_node);
  // nelements that will be stored in all receivers per node
  nelements_ram_receivers_node = nelements_ram_receivers_node / element_size;
  // nelements per receiver
  size_t nelements_ram_receiver = nelements_ram_receivers_node / my_info.nreceivers_node;
  my_info.nelements_mem_receiver = nelements_ram_receiver;

  return my_info;
}


/* Everyone should agree on the same seed at the beginning.
 * Update the arguments: mixer_seed, byte_hasher_seed */
void seed_agreement(MITM_MPI_data& my_info,
		    u64& mixer_seed,
		    u64& byte_hasher_seed)
{
  u64 seeds[2] = {0};
  // This is sender No. 0, assuming round robin distribution of ranks.
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


}

#endif
