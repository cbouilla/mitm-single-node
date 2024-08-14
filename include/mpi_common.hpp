#ifndef MITM_MPI_COMMON
#define MITM_MPI_COMMON

#include <mpi.h>
#include <err.h>

#include "common.hpp"
#include "engine_common.hpp"


namespace mitm {

enum tags {TAG_INTERCOMM, TAG_POINTS, TAG_SENDER_CALLHOME, TAG_RECEIVER_CALLHOME, TAG_ASSIGNMENT, TAG_SOLUTION};
enum role {CONTROLLER, SENDER, RECEIVER, UNDECIDED};
enum assignment {KEEP_GOING, NEW_VERSION};


class MpiCounters : public BaseCounters {
public:
    MpiCounters() {}
    MpiCounters(bool display_active) : BaseCounters(display_active) {}

	/* waiting times for MPI implementation */
	double send_wait = 0.;          // time sender processes wait for the receiver to be ready
	double recv_wait = 0.;          // time the receivers wait for data from the senders
	u64 bytes_sent = 0;

	/* display management */
	u64 bytes_sent_prev = 0;
	
	void reset()
	{
		n_points = 0;
		n_dp = 0;
		n_collisions = 0;
		send_wait = 0;
		recv_wait = 0;
		bytes_sent = 0;
	}

	void display()
	{
		if (not display_active)
			return;
		double now = wtime();
		double delta = now - last_display;
		if (delta >= min_delay) {
			u64 N = 1ull << pb_n;
			double dprate = (n_dp - n_dp_prev) / delta;
			char hdprate[8];
			human_format(dprate, hdprate);
			double frate = (n_points - n_points_prev) / delta;
			char hfrate[8];
			human_format(frate, hfrate);
			double nrate = (bytes_sent - bytes_sent_prev) / delta;
			char hnrate[8];
			human_format(nrate, hnrate);
			printf("\r#i = %" PRId64 " (%.02f*n/w). %s f/sec.  %s DP/sec.  #DP (this i / total) %.02f*w / %.02f*n.  #coll (this i / total) %.02f*w / %0.2f*n.  Total #f=2^%.02f.  S/R wait %.02fs / %.02fs.  S->R %sB/s",
			 		n_flush, (double) n_flush * w / N, 
			 		hfrate, hdprate, 
			 		(double) n_dp_i / w, (double) n_dp / N,
			 		(double) n_collisions_i / w, (double) n_collisions / N,
			 		std::log2(n_points),
			 		send_wait, recv_wait,
			 		hnrate);
			fflush(stdout);
			n_dp_prev = n_dp;
			n_points_prev = n_points;
			bytes_sent_prev = bytes_sent;
			last_display = wtime();
		}
	}
};  


class MpiParameters : public Parameters {
public:
	int recv_per_node = 1;
	int buffer_capacity = 1500;            // somewhat arbitrary
	double ping_delay = 0.1;

	MPI_Comm world_comm;
	MPI_Comm inter_comm;
	int role = UNDECIDED;                   /* enum role */
	int local_rank = 0;                     /* rank among the local group of the inter-communicator */
	int n_recv;
	int n_send;

	void setup(MPI_Comm comm)
	{
		world_comm = comm;
		int size, rank;
		MPI_Comm_size(world_comm, &size);
		MPI_Comm_rank(world_comm, &rank);
		verbose = 0;
		if (rank == 0) {
			role = CONTROLLER;
			verbose = 1;
		}
	
		/* create a subcommunicator inside each node */
		MPI_Comm node_comm;
		MPI_Comm_split_type(world_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
	
		/* determine the number of nodes (== #processes of node-rank 0) */
		int node_rank;
		int node_size;
		MPI_Comm_rank(node_comm, &node_rank);
		MPI_Comm_size(node_comm, &node_size);
		int is_rank0 = (node_rank == 0) ? 1 : 0;
		int n_nodes;
		MPI_Allreduce(&is_rank0, &n_nodes, 1, MPI_INT, MPI_SUM, world_comm);
		if (role == CONTROLLER) {
			printf("MPI: detected %d nodes\n", n_nodes);
			/* verification */
			if ((size % n_nodes) != 0)
				errx(1, "MPI: ERROR! The number of processes (%d) is not a multiple of the number of hosts (%d)", size, n_nodes);
		}
	
		/* check that there is always the same number of processes per node */
		int check_hi, check_lo;
		MPI_Allreduce(&node_size, &check_lo, 1, MPI_INT, MPI_MIN, world_comm);
		MPI_Allreduce(&node_size, &check_hi, 1, MPI_INT, MPI_MAX, world_comm);
		if (role == CONTROLLER) {
			if (check_lo != check_hi)
				printf("MPI: WARNING!!! process / node varies from %d to %d\n", check_lo, check_hi);
			else
				printf("MPI: each node runs %d processes\n", check_lo);
		}
		MPI_Comm_free(&node_comm);

#if 0 
		// NUMA stuff. Let's keep this for later
		/* split each node comm into a NUMA comm */
		MPI_Comm numa_comm;
		// this is according to the spec
		// MPI_Info info;
		// MPI_Info_create(&info);
		// MPI_Info_set(info, "mpi_hw_resource_type", "NUMANode");
		// MPI_Comm_split_type(node_comm, MPI_COMM_TYPE_HW_GUIDED, 0, info, &numa_comm);
	
		// this works for OpenMPI
		MPI_Comm_split_type(node_comm, OMPI_COMM_TYPE_NUMA, 0, MPI_INFO_NULL, &numa_comm);
	
		/* how many of us in the NUMA comm? */
		int numa_rank;
		int numa_size;
		MPI_Comm_size(numa_comm, &numa_size);
		MPI_Comm_rank(numa_comm, &numa_rank);
	
		/* count NUMA comms */
		int n_numa;
		int is_numarank0 = (numa_rank == 0) ? 1 : 0;
		MPI_Allreduce(&is_numarank0, &n_numa, 1, MPI_INT, MPI_SUM, world_comm);
	
		/* check size of NUMA comms */
		MPI_Allreduce(&numa_size, &check_lo, 1, MPI_INT, MPI_MIN, world_comm);
		MPI_Allreduce(&numa_size, &check_hi, 1, MPI_INT, MPI_MAX, world_comm);
	
		if (master) {
			if (check_lo != check_hi)
				printf("MPI: WARNING! #process / #NUMA zones varies from %d to %d\n", check_lo, check_hi);
			else
				printf("MPI: detected %d NUMA zones per node\n", size / check_lo / n_nodes);
		}
#endif	
		
		/* decide sender / receiver */
		if (node_rank >= node_size - recv_per_node)
			role = RECEIVER;
		else if (role == UNDECIDED)
			role = SENDER;
		/* safety check */
		assert(role != UNDECIDED);
		if (rank == 0 && role != CONTROLLER)
			errx(1, "something went wrong");
		/* count them */
		n_recv = (role == RECEIVER) ? 1 : 0;
		MPI_Allreduce(MPI_IN_PLACE, &n_recv, 1, MPI_INT, MPI_SUM, world_comm);
		n_send = size - n_recv - 1;           // exclude the CONTROLLER
		if (role == CONTROLLER) {
			printf("MPI: # sender   processes = %d\n", n_send);
			printf("MPI: # receiver processes = %d\n", n_recv);
			if (n_send == 0 || n_recv == 0)
				errx(1, "MPI: ERROR! At least one sender/receiver is required");
		}
	
		/* locate last receiver / sender */
		int last_rank[2];
		last_rank[0] = (role == SENDER) ? rank : 0;
		last_rank[1] = (role == RECEIVER) ? rank : 0;
		MPI_Allreduce(MPI_IN_PLACE, last_rank, 2, MPI_INT, MPI_MAX, world_comm);

		/* Create an intra-comm that separates senders and receivers */
		MPI_Comm tribe_comm;
		MPI_Comm_split(world_comm, role, 0, &tribe_comm);
		MPI_Comm_rank(tribe_comm, &local_rank);
		if (role != CONTROLLER) {
			/* creation of an inter-communicator between senders and receivers */
			int local_leader;     // in tribe_comm
			int remote_leader;    // in world_comm
			if (role == RECEIVER) {
				local_leader = n_recv - 1;
				remote_leader = last_rank[0];
			} else {
				local_leader = n_send - 1;
				remote_leader = last_rank[1];
			}
			MPI_Intercomm_create(tribe_comm, local_leader, world_comm, remote_leader, TAG_INTERCOMM, &inter_comm);
		}
		MPI_Comm_free(&tribe_comm);
	}
};

}
#endif