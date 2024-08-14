#ifndef MITM_MPI_RECEIVER
#define MITM_MPI_RECEIVER

#include <vector>
#include <mpi.h>

#include "mpi_common.hpp"
#include "engine_common.hpp"

namespace mitm {

/* Manage reception buffers for a collection of sender processes, with double-buffering */
class RecvBuffers {
	using Buffer = std::vector<u64>;

private:
	const MpiParameters &params;
	const size_t capacity;
	
	std::vector<Buffer> ready;                 // buffers containing points ready to be processed 
	std::vector<Buffer> incoming;              // buffers waiting for incoming data
	std::vector<MPI_Request> request;
	int n_active_senders;                      // # active senders


	/* initiate reception for a specific sender */
	void listen_sender(int i)
	{
		//printf("recv:listen_sender(%d)\n", i);
		assert(request[i] == MPI_REQUEST_NULL);
		incoming[i].resize(3 * capacity);
		MPI_Irecv(incoming[i].data(), 3 * capacity, MPI_UINT64_T, i, TAG_POINTS, params.inter_comm, &request[i]);
	}

public:
	RecvBuffers(const MpiParameters &params) : params(params), capacity(params.buffer_capacity)
	{
		int n = params.n_send;
		ready.resize(n);      /* does this even work? */
		incoming.resize(n);
		request.resize(n, MPI_REQUEST_NULL);
		for (int i = 0; i < n; i++) {
			ready[i].reserve(3 * capacity);
			incoming[i].reserve(3 * capacity);
		}
	}

	// workflow : (start (!complete wait)* complete shutdown)*

	/* initiate reception for all senders. Invoke at beginning, or after shutdown() */
	void start()
	{
		for (int i = 0; i < params.n_send; i++)
			listen_sender(i);
		n_active_senders = params.n_send;
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
		int n = params.n_send;
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

template<typename ConcreteProblem>
void receiver(const ConcreteProblem& Pb, const MpiParameters &params)
{
	RecvBuffers recvbuf(params);
    MpiCounters &ctr = Pb.ctr; 
    Dict<std::pair<u64, u64>> dict(params.nbytes_memory / params.recv_per_node);

    double last_ping = wtime();

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop	
		u64 i = msg[0];

		// receive and process data from senders
		recvbuf.start();
		for (;;) {
			if (recvbuf.complete())
				break;                      // all senders are done
			auto ready = recvbuf.wait(ctr);
			// process incoming buffers of distinguished points
			for (auto it = ready.begin(); it != ready.end(); it++) {
				auto & buffer = *it;
				for (size_t k = 0; k < buffer->size(); k += 3) {
					u64 start = buffer->at(k);
					u64 end = buffer->at(k + 1);
					u64 len = buffer->at(k + 2);
					assert((start & Pb.mask) == start);
					assert((end & Pb.mask) == end);
					auto solution = process_distinguished_point(Pb, ctr, dict, i, start, end, len);
					if (solution) {          // call home !
						auto [i, x0, x1] = *solution;
						u64 golden[3] = {i, x0, x1};
						MPI_Send(golden, 3, MPI_UINT64_T, 0, TAG_SOLUTION, params.world_comm);
					}
				}
			}
			
			/* call home? */
			if (wtime() - last_ping >= 1) {
            	last_ping = wtime();
            	u64 stats[3] = {ctr.n_points, ctr.n_collisions, (u64) (ctr.recv_wait * 1e6)};
            	// todo: send stats:  [#bytes received?]
            	MPI_Send(stats, 3, MPI_UINT64_T, 0, TAG_RECEIVER_CALLHOME, params.world_comm);
            	ctr.reset();
            }
		}
		dict.flush();
	}
}

}
#endif
