#ifndef MITM_MPI_COMMON
#define MITM_MPI_COMMON

#include <mpi.h>
#include <err.h>

#include "../common.hpp"

namespace mitm {

// enum tags {TAG_INTERCOMM, TAG_POINTS, TAG_SENDER_CALLHOME, TAG_RECEIVER_CALLHOME, TAG_ASSIGNMENT, TAG_SOLUTION};
// enum role {CONTROLLER, SENDER, RECEIVER, UNDECIDED};
// enum assignment {KEEP_GOING, NEW_VERSION};


class MpiParameters : public Parameters {
public:
	double ping_delay = 0.1;
	MPI_Comm comm;
	int mpi_rank, mpi_size;

	size_t buffer_capacity = 4096;               // somewhat arbitrary
	int n_recv_buffers = 32;                     // somewhat arbitrary
	int n_slack_buffers = 256;                   // somewhat arbitrary
    double send_buffers_ratio = 1.25;            // allocated send buffers = this * strictly required
    int min_backoff = 1000;                      // µs
    int max_backoff = 100000;                    // µs
    int max_current_send = 16;                   // somewhat arbitrary

	MpiParameters(MPI_Comm comm) : comm(comm)
	{
		MPI_Comm_rank(comm, &mpi_rank);
		MPI_Comm_size(comm, &mpi_size);
		verbose = (mpi_rank == 0);
	}
};

void BCast_result(const MpiParameters &params, vector<pair<u64,u64>> &result)
{
	// deal with the results
    vector<int> recvcounts(params.mpi_size);
    recvcounts[params.mpi_rank] = 2 * result.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, params.comm);

    vector<int> displs(params.mpi_size);
    int acc = 0;
    for (int i = 0; i < params.mpi_size; i++) {
        displs[i] = acc;
        acc += recvcounts[i];
    }

    vector<u64> tmp(acc);
    for (size_t i = 0; i < result.size(); i++) {
        int offset = displs[params.mpi_rank] + 2*i;
        auto [x0, x1] = result[i];
        tmp[offset] = x0;
        tmp[offset + 1] = x1;
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT, tmp.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, params.comm);
    result.clear();
    for (int i = 0; i < acc; i += 2)
        result.push_back(pair(tmp[i], tmp[i+1]));
}

}
#endif