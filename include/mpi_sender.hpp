#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER
#include <err.h>
#include <vector>

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "omp_engine.hpp"
#include "mpi_common.hpp"

namespace mitm {



/* Manages send buffers for a collection of receiver processes, with double-buffering */
class SendBuffers {
	using Buffer = std::vector<u64>;

private:
	const MpiParameters &params;
	const size_t capacity;
	
	std::vector<Buffer> ready;
	std::vector<Buffer> outgoing;
	std::vector<MPI_Request> request;   /* for the OUTGOING buffers */

	/* wait until the i-th passive buffer has been fully sent */
	void wait_send(int i)
	{
		MPI_Wait(&request[i], MPI_STATUS_IGNORE);
		outgoing[i].clear();
	}

	/* initiate transmission of the i-th OUTGOING buffer */
	void start_send(int i)
	{
		MPI_Issend(outgoing[i].data(), 3 * outgoing[i].size(), MPI_UINT64_T, i, TAG_POINTS, params.inter_comm, &request[i]);
	}

public:
	SendBuffers(const MpiParameters &params) : params(params), capacity(params.buffer_capacity)
	{
		int n = params.n_recv;
		ready.resize(n);      /* does this even work? */
		outgoing.resize(n);
		request.resize(n, MPI_REQUEST_NULL);
		for (int i = 0; i < n; i++) {
			ready[i].reserve(3*capacity);
			outgoing[i].reserve(3*capacity);
		}
	}

	/* add a new item to the send buffer. Send if necessary */
	void push(u64 start, u64 end, u64 len, int rank)
	{
		if (ready[rank].size() == 3*capacity) {
			/* ready buffer is full. */
			wait_send(rank);      // finish sending the outgoing buffer
			std::swap(ready[rank], outgoing[rank]);
			start_send(rank);     // start sending
		}
		ready[rank].push_back(start);
		ready[rank].push_back(end);
		ready[rank].push_back(len);
	}

	/* send and empty all buffers, even if they are incomplete */
	void flush()
	{
		int n = params.n_recv;
		// finish sending all the outgoing buffers
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// send all the (incomplete) ready buffers
		for (int i = 0; i < n; i++) {
			std::swap(ready[i], outgoing[i]);
			start_send(i);
		}
		MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

		// finally tell all receivers that we are done
		for (int i = 0; i < n; i++)
			MPI_Send(NULL, 0, MPI_INT, i, TAG_NO_MORE_DATA, params.inter_comm);
	}
};


template<typename ConcreteProblem>
void sender(const ConcreteProblem& Pb, const MpiParameters &params)
{
	SendBuffers sendbuf(params);
    Counters &ctr = Pb.ctr; 

    u64 steps = 0;
    double last_ping = wtime();

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop

		u64 i = msg[0];
		for (u64 start = msg[1] + params.local_rank;; start += params.n_send) {
			steps += 1;

			/* start a new chain from a fresh "random" starting point */
            auto dp = generate_dist_point(Pb, i, params, start & Pb.mask);
            if (not dp) {
                ctr.dp_failure();
                continue;       /* bad chain start */
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);

            u64 hash = (end * 0xdeadbeef) % 0x7fffffff;
            int target_recv = ((int) hash) % params.n_recv;
            sendbuf.push(start, end, len, target_recv);
		
            if ((steps >= 1000) && (wtime() - last_ping >= 1)) {
				/* call home */
				last_ping = wtime();
            	u64 buffer[1] = {ctr.n_dp_i};
            	MPI_Send(buffer, 1, MPI_UINT64_T, 0, TAG_SENDER_CALLHOME, params.world_comm);
            	int assignment;
            	MPI_Recv(&assignment, 1, MPI_INT, 0, TAG_ASSIGNMENT, params.world_comm, MPI_STATUS_IGNORE);
            	if (assignment == KEEP_GOING)
            		continue;
            	sendbuf.flush();   
            	if (assignment == NEW_VERSION)        /* new broadcast */
            		break;
            }
		}

		ctr.flush_dict();
	}
}

}
#endif