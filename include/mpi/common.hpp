#ifndef MITM_MPI_COMMON
#define MITM_MPI_COMMON

#include <mpi.h>
#include <err.h>

#include "../common.hpp"

namespace mitm {

enum tags {TAG_INTERCOMM, TAG_POINTS, TAG_SENDER_CALLHOME, TAG_RECEIVER_CALLHOME, TAG_ASSIGNMENT, TAG_SOLUTION};
enum role {CONTROLLER, SENDER, RECEIVER, UNDECIDED};
enum assignment {KEEP_GOING, NEW_VERSION};


class MpiParameters : public Parameters {
public:
	int recv_per_node = 1;
	int buffer_capacity = 1500;            // somewhat arbitrary
	double ping_delay = 0.1;

	MPI_Comm world_comm;
	MPI_Comm inter_comm;
	MPI_Comm local_comm;                    /* just our side of the intercomm */
	int role = UNDECIDED;                   /* enum role */
	int rank, size;                         // for the global communicator
	int local_rank, local_size;             /* rank among the local group of the inter-communicator */
	int n_send;
	int n_nodes;

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
		MPI_Comm_split(world_comm, role, 0, &local_comm);
		MPI_Comm_rank(local_comm, &local_rank);
		MPI_Comm_size(local_comm, &local_size);
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
			MPI_Intercomm_create(local_comm, local_leader, world_comm, remote_leader, TAG_INTERCOMM, &inter_comm);
		}
	}
};


/* Manages send buffers for a collection of receiver processes, with double-buffering */
class SendBuffers {
public:
	using Buffer = vector<u64>;
	double waiting_time = 0;
	u64 bytes_sent = 0;

private:
	MPI_Comm inter_comm;
	size_t capacity;
	int tag;
	int n;
	
	vector<Buffer> ready;
	vector<Buffer> outgoing;
	vector<MPI_Request> request;   /* for the OUTGOING buffers */

	/* initiate transmission of the i-th OUTGOING buffer */
	void start_send(int i)
	{
		if (outgoing[i].size() == 0)  // do NOT send empty buffers. These are interpreted as "I am done"
			return;
		MPI_Isend(outgoing[i].data(), outgoing[i].size(), MPI_UINT64_T, i, tag, inter_comm, &request[i]);
		bytes_sent += outgoing[i].size() * sizeof(u64);
	}

	void switch_when_full(int rank)
	{
		if (ready[rank].size() == capacity) {
			// ready buffer is full: finish sending the outgoing buffer
			double start = wtime();
			MPI_Wait(&request[rank], MPI_STATUS_IGNORE);
			waiting_time += wtime() - start;
			outgoing[rank].clear();
			// swap and start sending
			std::swap(ready[rank], outgoing[rank]);
			start_send(rank);
		}
	}

public:
	SendBuffers(MPI_Comm inter_comm, int tag, size_t capacity) : inter_comm(inter_comm), capacity(capacity), tag(tag)
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
	void push(u64 x, int rank)
	{
		switch_when_full(rank);
		ready[rank].push_back(x);
	}

	void push2(u64 x, u64 y, int rank)
	{
		switch_when_full(rank);
		ready[rank].push_back(x);
		ready[rank].push_back(y);
	}

	void push3(u64 x, u64 y, u64 z, int rank)
	{
		switch_when_full(rank);
		ready[rank].push_back(x);
		ready[rank].push_back(y);
		ready[rank].push_back(z);
	}

	/* send and empty all buffers, even if they are incomplete */
	void flush()
	{
		// finish sending all the outgoing buffers
		double start = wtime();
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// send all the (incomplete) ready buffers
		for (int i = 0; i < n; i++) {
			std::swap(ready[i], outgoing[i]);
			start_send(i);
		}
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// finally tell all receivers that we are done
		for (int i = 0; i < n; i++)
			MPI_Send(NULL, 0, MPI_UINT64_T, i, tag, inter_comm);
		waiting_time += wtime() - start;
	}
};


/* Manage reception buffers for a collection of sender processes, with double-buffering */
class RecvBuffers {
public:
	using Buffer = vector<u64>;
	double waiting_time = 0;
	u64 bytes_sent = 0;

private:
	MPI_Comm inter_comm;
	const size_t capacity;
	int n;
	int tag;

	vector<Buffer> ready;                 // buffers containing points ready to be processed 
	vector<Buffer> incoming;              // buffers waiting for incoming data
	vector<MPI_Request> request;

	/* initiate reception for a specific sender */
	void listen_sender(int i)
	{
		incoming[i].resize(capacity);
		MPI_Irecv(incoming[i].data(), capacity, MPI_UINT64_T, i, tag, inter_comm, &request[i]);
	}

public:
	int n_active_senders;                      // # active senders
	RecvBuffers(MPI_Comm inter_comm, int tag, size_t capacity) : inter_comm(inter_comm), capacity(capacity), tag(tag)
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

	~RecvBuffers()
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
	vector<Buffer *> wait()
	{
		assert(n_active_senders > 0);
		vector<Buffer *> result;
		int n_done;
		vector<int> rank_done(n);
		vector<MPI_Status> statuses(n);
		double start = wtime();
		MPI_Waitsome(n, request.data(), &n_done, rank_done.data(), statuses.data());
		waiting_time += wtime() - start;
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

void BCast_result(MpiParameters &params, vector<pair<u64,u64>> &result)
{
	// deal with the results
    vector<int> recvcounts(params.size);
    recvcounts[params.rank] = 2 * result.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, params.world_comm);

    vector<int> displs(params.size);
    int acc = 0;
    for (int i = 0; i < params.size; i++) {
        displs[i] = acc;
        acc += recvcounts[i];
    }

    vector<u64> tmp(acc);
    for (size_t i = 0; i < result.size(); i++) {
        int offset = displs[params.rank] + 2*i;
        auto [x0, x1] = result[i];
        tmp[offset] = x0;
        tmp[offset + 1] = x1;
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT, tmp.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, params.world_comm);
    result.clear();
    for (int i = 0; i < acc; i += 2)
        result.push_back(pair(tmp[i], tmp[i+1]));
}

}
#endif