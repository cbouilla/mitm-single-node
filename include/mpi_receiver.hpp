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
	std::vector<MPI_Request> request_points;   // for the incomming buffers
	std::vector<MPI_Request> request_nomore;
	int n_active_senders;                      // # active senders


	/* initiate reception for a specific sender */
	void listen_sender(int i)
	{
		assert(request_points[i] == MPI_REQUEST_NULL);
		incoming[i].resize(3 * capacity);
		MPI_Irecv(incoming[i].data(), 3 * capacity, MPI_UINT64_T, i, TAG_POINTS, params.inter_comm, &request_points[i]);
	}

public:
	RecvBuffers(const MpiParameters &params) : params(params), capacity(params.buffer_capacity)
	{
		int n = params.n_send;
		ready.resize(n);      /* does this even work? */
		incoming.resize(n);
		request_points.resize(n, MPI_REQUEST_NULL);
		request_nomore.resize(n, MPI_REQUEST_NULL);
		for (int i = 0; i < n; i++) {
			ready[i].reserve(3 * capacity);
			incoming[i].reserve(3 * capacity);
		}
	}

	// workflow : (start (!complete wait)* complete shutdown)*

	/* initiate reception for all senders. Invoke at beginning, or after shutdown() */
	void start()
	{
		for (int i = 0; i < params.n_send; i++) {
			listen_sender(i);
			assert(request_nomore[i] == MPI_REQUEST_NULL);
			MPI_Irecv(NULL, 0, MPI_INT, i, TAG_NO_MORE_DATA, params.inter_comm, &request_nomore[i]);
		}
	}

	/* 
	 * invoke when complete() returns true to clean everything up.
	 * All requests are inactive when this terminates
	 */
	void shutdown()
	{
		assert(n_active_senders == 0);
		// The request_points[...] have all been cancelled
		for (int i = 0; i < params.n_send; i++)
			assert(request_nomore[i] == MPI_REQUEST_NULL);
		
		// clear the requests
		std::vector<MPI_Status> statuses(params.n_send);
		MPI_Waitall(params.n_send, request_points.data(), statuses.data());
		// verify that cancellation was succesful
		for (int i = 0; i < params.n_send; i++) {
			int cancelled;
			MPI_Test_cancelled(&statuses[i], &cancelled);
			assert(cancelled);
		}
	}

	/* return true when all senders are done; update n_active_senders */
	bool complete()
	{
		int n = params.n_send;
		int n_done;
		std::vector<int> rank_done(n);
		MPI_Testsome(n, request_nomore.data(), &n_done, rank_done.data(), MPI_STATUSES_IGNORE);
		if (n_done == MPI_UNDEFINED) {
			assert(n_active_senders == 0);
			return 1;     // no active requests ---> no sender is going to send anything anymore
		}
		for (int i = 0; i < n_done; i++) {
			MPI_Wait(&request_nomore[i], MPI_STATUS_IGNORE);
			MPI_Cancel(&request_points[i]);
		}
		n_active_senders -= n_done;
		return (n_active_senders == 0);
	}

	/* 
	 * Wait until some data arrives. Returns the buffers that have arrived.
	 * Only call this when complete() returned false (otherwise, this will wait forever)
	 * This may destroy the content of all "ready" buffers, so that they have to be processed first
	 */
	std::vector<Buffer *> wait()
	{
		assert(n_active_senders > 0);
		int n = params.n_send;
		std::vector<Buffer *> result;
		int n_done;
		std::vector<int> rank_done(n);
		std::vector<MPI_Status> statuses(n);
		MPI_Waitsome(n, request_points.data(), &n_done, rank_done.data(), statuses.data());

		for (int i = 0; i < n_done; i++) {
			int j = rank_done[i];
			std::swap(incoming[j], ready[j]);
			int count;
			MPI_Get_count(&statuses[j], MPI_UINT64_T, &count);
			ready[j].resize(count);         // matching message size
			result.push_back(&ready[j]);
			listen_sender(j);
		}
		return result;
	}
};

template<typename ConcreteProblem>
void receiver(const ConcreteProblem& Pb, const MpiParameters &params)
{
	RecvBuffers recvbuf(params);
    Counters &ctr = Pb.ctr; 
    Dict<std::pair<u64, u64>> dict(params.nbytes_memory / params.recv_per_node);

    u64 steps = 0;
    double last_ping = wtime();

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop
		
		u64 i = msg[0];
		recvbuf.start();
		for (;;) {
			if (recvbuf.complete())
				break;                      // senders are done
			// wait for senders
			auto ready = recvbuf.wait();
			// process incoming buffers of distinguished points
			for (auto it = ready.begin(); it != ready.end(); it++) {
				auto & buffer = *it;
				for (int i = 0; i < buffer->size(); i += 3) {
					u64 start = buffer->at(i);
					u64 end = buffer->at(i + 1);
					u64 len = buffer->at(i + 2);
					auto solution = process_distinguished_point(Pb, ctr, dict, i, start, end, len);
					if (solution) {          // call home !
						u64 golden[3] = *solution;
						MPI_Send(golden, 3, MPI_UINT64_T, 0, TAG_SOLUTION, params.world_comm);
					}
				}
			}
			
			/* call home? */
			if (wtime() - last_ping >= 1) {
            	last_ping = wtime();
            	// todo: send stats?
            	MPI_Send(NULL, 0, MPI_INT, 0, TAG_RECEIVER_CALLHOME, params.world_comm);
            }
		}
		recvbuf.shutdown();
	}
}

}
#endif
