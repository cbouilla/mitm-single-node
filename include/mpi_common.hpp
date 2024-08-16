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
	int rank, size;                         // for the global communicator
	int local_rank = 0;                     /* rank among the local group of the inter-communicator */
	int n_recv;
	int n_send;

	void setup(MPI_Comm comm)
	{
		setup(comm, 1);
	}

	void setup(MPI_Comm comm, bool controller)
	{
		world_comm = comm;
		MPI_Comm_size(world_comm, &size);
		MPI_Comm_rank(world_comm, &rank);
		verbose = (rank == 0);
		if (controller && rank == 0)
			role = CONTROLLER;
	
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
		if (verbose) {
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
		/* count them */
		n_recv = (role == RECEIVER) ? 1 : 0;
		MPI_Allreduce(MPI_IN_PLACE, &n_recv, 1, MPI_INT, MPI_SUM, world_comm);
		n_send = (role == SENDER) ? 1 : 0;
		MPI_Allreduce(MPI_IN_PLACE, &n_send, 1, MPI_INT, MPI_SUM, world_comm);
		if (verbose) {
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


/* Manages send buffers for a collection of receiver processes, with double-buffering */
class BaseSendBuffers {
public:
	using Buffer = std::vector<u64>;

protected:
	MPI_Comm inter_comm;
	size_t capacity;
	int tag;
	int n;
	
	std::vector<Buffer> ready;
	std::vector<Buffer> outgoing;
	std::vector<MPI_Request> request;   /* for the OUTGOING buffers */

	/* wait until the i-th passive buffer has been fully sent */
	void wait_send(int i, MpiCounters &ctr)
	{
		
	}

	/* initiate transmission of the i-th OUTGOING buffer */
	void start_send(int i, MpiCounters &ctr)
	{
		if (outgoing[i].size() == 0)  // do NOT send empty buffers. These are interpreted as "I am done"
			return;
		MPI_Issend(outgoing[i].data(), outgoing[i].size(), MPI_UINT64_T, i, tag, inter_comm, &request[i]);
		ctr.bytes_sent += outgoing[i].size() * sizeof(u64);
	}

	void switch_when_full(int rank, MpiCounters &ctr)
	{
		if (ready[rank].size() == capacity) {
			// ready buffer is full: finish sending the outgoing buffer
			double start = wtime();
			MPI_Wait(&request[rank], MPI_STATUS_IGNORE);
			ctr.send_wait += wtime() - start;
			outgoing[rank].clear();
			// swap and start sending
			std::swap(ready[rank], outgoing[rank]);
			start_send(rank, ctr);
		}
	}

public:
	BaseSendBuffers(MPI_Comm inter_comm, int tag, size_t capacity) : inter_comm(inter_comm), capacity(capacity), tag(tag)
	{
		MPI_Comm_remote_size(inter_comm, &n);
		ready.resize(n);
		outgoing.resize(n);
		request.resize(n, MPI_REQUEST_NULL);
		for (int i = 0; i < n; i++) {
			ready[i].reserve(capacity);
			outgoing[i].reserve(capacity);
		}
	}

	/* add a new item to the send buffer. Send if necessary */
	void push(u64 x, int rank, MpiCounters &ctr)
	{
		switch_when_full(rank, ctr);
		ready[rank].push_back(x);
	}

	/* send and empty all buffers, even if they are incomplete */
	void flush(MpiCounters &ctr)
	{
		// finish sending all the outgoing buffers
		double start = wtime();
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// send all the (incomplete) ready buffers
		for (int i = 0; i < n; i++) {
			std::swap(ready[i], outgoing[i]);
			start_send(i, ctr);
		}
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// finally tell all receivers that we are done
		for (int i = 0; i < n; i++)
			MPI_Send(NULL, 0, MPI_UINT64_T, i, tag, inter_comm);
		ctr.send_wait += wtime() - start;
	}
};


/* Manage reception buffers for a collection of sender processes, with double-buffering */
class BaseRecvBuffers {
	using Buffer = std::vector<u64>;

private:
	MPI_Comm inter_comm;
	const size_t capacity;
	int n;
	int tag;

	std::vector<Buffer> ready;                 // buffers containing points ready to be processed 
	std::vector<Buffer> incoming;              // buffers waiting for incoming data
	std::vector<MPI_Request> request;
	int n_active_senders;                      // # active senders


	/* initiate reception for a specific sender */
	void listen_sender(int i)
	{
		incoming[i].resize(capacity);
		MPI_Irecv(incoming[i].data(), capacity, MPI_UINT64_T, i, tag, inter_comm, &request[i]);
	}

public:
	BaseRecvBuffers(MPI_Comm inter_comm, int tag, size_t capacity) : inter_comm(inter_comm), capacity(capacity), tag(tag)
	{
		MPI_Comm_remote_size(inter_comm, &n);
		ready.resize(n);
		incoming.resize(n);
		request.resize(n, MPI_REQUEST_NULL);
		for (int i = 0; i < n; i++) {
			ready[i].reserve(capacity);
			incoming[i].reserve(capacity);
			listen_sender(i);
		}
		n_active_senders = n;
	}

	~BaseRecvBuffers()
	{
		assert (n_active_senders == 0);
	}

	/* return true when all senders are done */
	bool complete()
	{
		return (n_active_senders == 0);
	}

	/* 
	 * Wait until some data arrives. Returns the buffers that have arrived.
	 * Only call this when complete() returned false (otherwise, this will wait forever)
	 * This may destroy the content of all "ready" buffers, so that they have to be processed first
	 */
	std::vector<Buffer *> wait(MpiCounters &ctr)
	{
		assert(n_active_senders > 0);
		std::vector<Buffer *> result;
		int n_done;
		std::vector<int> rank_done(n);
		std::vector<MPI_Status> statuses(n);
		double start = wtime();
		MPI_Waitsome(n, request.data(), &n_done, rank_done.data(), statuses.data());
		ctr.recv_wait += wtime() - start;
		assert(n_done != MPI_UNDEFINED);

		for (int i = 0; i < n_done; i++) {
			int j = rank_done[i];
			std::swap(incoming[j], ready[j]);
			int count;
			MPI_Get_count(&statuses[i], MPI_UINT64_T, &count);
			if (count == 0) {
				n_active_senders -= 1;
			} else {
				ready[j].resize(count);         // matching message size
				result.push_back(&ready[j]);
				listen_sender(j);
			}
		}
		return result;
	}
};

}
#endif