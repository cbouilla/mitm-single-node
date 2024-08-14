#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER
#include <err.h>
#include <vector>

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
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
	void wait_send(int i, MpiCounters &ctr)
	{
		double start = wtime();
		MPI_Wait(&request[i], MPI_STATUS_IGNORE);
		ctr.send_wait += wtime() - start;
		outgoing[i].clear();
	}

	/* initiate transmission of the i-th OUTGOING buffer */
	void start_send(int i, MpiCounters &ctr)
	{
		if (outgoing[i].size() == 0)  // do NOT send empty buffers. These are interpreted as "I am done"
			return;
		MPI_Issend(outgoing[i].data(), outgoing[i].size(), MPI_UINT64_T, i, TAG_POINTS, params.inter_comm, &request[i]);
		ctr.bytes_sent += outgoing[i].size() * sizeof(u64);
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
	void push(u64 start, u64 end, u64 len, int rank, MpiCounters &ctr)
	{
		if (ready[rank].size() == 3 * capacity) {
			/* ready buffer is full. */
			wait_send(rank, ctr);      // finish sending the outgoing buffer
			std::swap(ready[rank], outgoing[rank]);
			start_send(rank, ctr);     // start sending
		}
		ready[rank].push_back(start);
		ready[rank].push_back(end);
		ready[rank].push_back(len);
	}

	/* send and empty all buffers, even if they are incomplete */
	void flush(MpiCounters &ctr)
	{
		int n = params.n_recv;
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
			MPI_Send(NULL, 0, MPI_UINT64_T, i, TAG_POINTS, params.inter_comm);
		ctr.send_wait += wtime() - start;
	}
};


template<typename ConcreteProblem>
void sender(const ConcreteProblem& Pb, const MpiParameters &params)
{
	SendBuffers sendbuf(params);
    MpiCounters &ctr = Pb.ctr; 

    u64 steps = 0;
    double last_ping = wtime();

	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop

		u64 i = msg[0];
		for (u64 j = msg[1] + 3*params.local_rank;; j += 3*params.n_send) {   // add an odd number to avoid problems mod 2^n...
			steps += 1;

			/* call home? */
            if ((steps % 10000 == 0) && (wtime() - last_ping >= params.ping_delay)) {
				steps = 0;
				last_ping = wtime();
            	u64 stats[4] = {ctr.n_points, ctr.n_dp, (u64) (ctr.send_wait * 1e6), ctr.bytes_sent};
            	MPI_Send(stats, 4, MPI_UINT64_T, 0, TAG_SENDER_CALLHOME, params.world_comm);
            	ctr.reset();

            	int assignment;
            	MPI_Recv(&assignment, 1, MPI_INT, 0, TAG_ASSIGNMENT, params.world_comm, MPI_STATUS_IGNORE);
            	if (assignment == NEW_VERSION) {        /* new broadcast */
            	   	sendbuf.flush(ctr);   
            		break;
            	}
            }

			/* start a new chain from a fresh "random" starting point */
			u64 start = j & Pb.mask;
            auto dp = generate_dist_point(Pb, i, params, start);
            if (not dp) {
                ctr.dp_failure();
                continue;       /* bad chain start */
            }

            auto [end, len] = *dp;
            ctr.found_distinguished_point(len);

            u64 hash = (end * 0xdeadbeef) % 0x7fffffff;
            int target_recv = ((int) hash) % params.n_recv;
            sendbuf.push(start, end, len, target_recv, ctr);
		}
	}
}

}
#endif