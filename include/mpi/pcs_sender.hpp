#ifndef MITM_MPI_SENDER
#define MITM_MPI_SENDER
#include <err.h>
#include <vector>

#include <mpi.h>

#include "common.hpp"
#include "engine_common.hpp"
#include "mpi/common.hpp"

namespace mitm {

// not DRY wrt sequential/pcs_engine.hpp
static void start_chain(const Parameters &params, u64 out_mask, u64 root_seed, u64 &j, u64 x[], u64 len[], u64 seed[], u64 jinc, int k)
{
    u64 start;
    for (;;) {
        j += jinc;
        start = (root_seed + j * params.multiplier) & out_mask;
        if (not is_distinguished_point(start, params.threshold))  // refuse to start from a DP
            break;
    }
    x[k] = start;
    len[k] = 0;
    seed[k] = j;    
}


template<class ProblemWrapper>
void sender(ProblemWrapper& wrapper, const MpiParameters &params)
{
    int jbits = std::log2(10 * params.w) + 8;
    u64 jmask = make_mask(jbits);
	for (;;) {
		/* get data from controller */
		u64 msg[3];   // i, root_seed, stop?
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.world_comm);
		if (msg[2] != 0)
			return;      // controller tells us to stop

    	u64 n_dp = 0;    // #DP found since last report
    	wrapper.n_eval = 0;
		SendBuffers sendbuf(params.inter_comm, TAG_POINTS, 3 * params.buffer_capacity);
    	double last_ping = wtime();
		u64 i = msg[0];
		u64 root_seed = msg[1];

		/* current state of the chains */
		constexpr int vlen = ProblemWrapper::vlen;
    	u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
    	u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
    	u64 len[vlen], seed[vlen];
		u64 j = params.local_rank;

		/* infinite loop to generate DPs */
        for (int k = 0; k < vlen; k++)
            start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, params.n_send, k);
		assert((j & jmask) == j);

		for (;;) {
			/* call home? */
            if ((n_dp % 10000 == 9999) && (wtime() - last_ping >= params.ping_delay)) {
				last_ping = wtime();
            	MPI_Send(&n_dp, 1, MPI_UINT64_T, 0, TAG_SENDER_CALLHOME, params.world_comm);
				n_dp = 0;

            	int assignment;
            	MPI_Recv(&assignment, 1, MPI_INT, 0, TAG_ASSIGNMENT, params.world_comm, MPI_STATUS_IGNORE);
            	if (assignment == NEW_VERSION) {        /* new broadcast */
            	   	sendbuf.flush();   
            		break;
            	}
            }

			/* advance all the chains */
        	wrapper.vmixf(i, x, y);

			/* test for distinguished points */ 
			for (int k = 0; k < vlen; k++) {
			    len[k] += 1;
			    x[k] = y[k];
			    bool dp = is_distinguished_point(x[k], params.threshold);
			    bool failure = (len[k] == params.dp_max_it);
			    if (dp) {
					n_dp += 1;
					int target_recv = (int) (x[k] % params.n_recv);
					sendbuf.push3(seed[k], x[k] / params.n_recv, len[k], target_recv);			        
			    }
			    if (dp || failure) {
			        start_chain(params, wrapper.out_mask, root_seed, j, x, len, seed, params.n_send, k);
			        assert((j & jmask) == j);
			    }
			}
		}

		// now is a good time to collect stats
		//             #f send,   
		u64 iavg[7] = {wrapper.n_eval, 0, 0, 0, 0, 0, 0};
		MPI_Reduce(iavg, NULL, 7, MPI_UINT64_T, MPI_SUM, 0, params.world_comm);
		//                send wait             recv wait
		double dmin[2] = {sendbuf.waiting_time, HUGE_VAL};
		double dmax[2] = {sendbuf.waiting_time, 0};
		double davg[2] = {sendbuf.waiting_time, 0};
		MPI_Reduce(dmin, NULL, 2, MPI_DOUBLE, MPI_MIN, 0, params.world_comm);
		MPI_Reduce(dmax, NULL, 2, MPI_DOUBLE, MPI_MAX, 0, params.world_comm);
		MPI_Reduce(davg, NULL, 2, MPI_DOUBLE, MPI_SUM, 0, params.world_comm);

		vector<u8> hll(0x10000);
		MPI_Reduce(hll.data(), NULL, 0x10000, MPI_UINT8_T, MPI_MAX, 0, params.world_comm);
	}
}

}
#endif