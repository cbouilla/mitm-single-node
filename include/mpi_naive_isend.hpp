#ifndef MITM_NAIVE_MPI_ISEND
#define MITM_NAIVE_MPI_ISEND

#include <optional>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cmath>

#include "common.hpp"
#include "AbstractCollisionProblem.hpp"
#include "mpi_common.hpp"

#include <mpi.h>

/*
 * naive MITM w/ distributed dictionnary.  round-based Alltoallv version.
 */

namespace mitm {

#define EXPENSIVE_F

template <class AbstractProblem>
std::vector<std::pair<u64, u64>> naive_mpi_claw_search_isend(AbstractProblem &Pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    double start = wtime();
    MpiCounters ctr;
    u64 N = 1ull << Pb.n;
    std::vector<std::pair<u64, u64>> result;
    CompactDict dict((params.role == RECEIVER) ? (1.25 * N) / params.n_recv : 0);

    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        if (params.verbose)
            printf("Starting phase %d\n", phase);

        double phase_start = wtime();
        ctr.reset();
        
        double wait;

        if (params.role == SENDER) {
            BaseSendBuffers sendbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            u64 lo = params.local_rank * N / params.n_send;
            u64 hi = (params.local_rank + 1) * N / params.n_send;
            for (u64 x = lo; x < hi; x++) {
                u64 z = (phase == 0) ? Pb.f(x) : Pb.g(x);
                u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
                int target = ((int) hash) % params.n_recv;
                #ifdef EXPENSIVE_F
                sendbuf.push2(x, z, target, ctr);
                #else
                sendbuf.push(x, target, ctr);
                #endif
            }
            sendbuf.flush(ctr);

            /* aggregate stats over all senders */
            wait = ctr.send_wait;
        }

        if (params.role == RECEIVER) {
            BaseRecvBuffers recvbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            u64 keys[3 * Pb.n];
            while (not recvbuf.complete()) {
                auto ready_buffers = recvbuf.wait(ctr);
                for (auto it = ready_buffers.begin(); it != ready_buffers.end(); it++) {
                    auto * buffer = *it;
                    // printf("got buffer! phase=%d, size=%zd\n", phase, buffer->size());
                    for (auto jt = buffer->begin(); jt != buffer->end(); jt++) {
                        u64 x = *jt;
                        #ifdef EXPENSIVE_F
                        jt++;
                        u64 z = *jt;
                        #else
                        u64 z = (phase == 0) ? Pb.f(x) : Pb.g(x);
                        #endif

                        switch (phase) {
                        case 0:
                            dict.insert(z, x);
                            break;
                        case 1: {
                            // probe dict
                            int nkeys = dict.probe(z, keys);
                            for (int k = 0; k < nkeys; k++) {
                                u64 y = keys[k];
                                if (z != Pb.f(y)) {
                                    ctr.collision_failure();
                                    continue;    // false positive from truncation in the hash table
                                }
                                ctr.found_collision();
                                if (Pb.is_good_pair(y, x)) {
                                    // printf("\nfound golden collision !!!\n");
                                    result.push_back(std::pair(y, x));
                                }
                            }
                        }
                        }
                    }
                }
            }
            wait = ctr.recv_wait;
        } // RECEIVER
        
        // timing
        MPI_Barrier(MPI_COMM_WORLD);
        double delta = wtime() - phase_start;
        double wait_min = wait;
        double wait_max = wait;
        double wait_avg = wait;
        MPI_Allreduce(MPI_IN_PLACE, &wait_min, 1, MPI_DOUBLE, MPI_MIN, params.local_comm);
        MPI_Allreduce(MPI_IN_PLACE, &wait_max, 1, MPI_DOUBLE, MPI_MAX, params.local_comm);
        MPI_Allreduce(MPI_IN_PLACE, &wait_avg, 1, MPI_DOUBLE, MPI_SUM, params.local_comm);
        wait_avg /= params.local_size;
        double wait_std = (wait - wait_avg) * (wait - wait_avg);
        MPI_Allreduce(MPI_IN_PLACE, &wait_std, 1, MPI_DOUBLE, MPI_SUM, params.local_comm);            
        wait_std = std::sqrt(wait_std);
        if (params.local_rank == 0) {
            printf("phase %d %s, wait min %.2fs max %.2fs avg %.2fs (%.1f%%) std %.2fs.\n",
                phase, (params.role == SENDER) ? "sender" : "receiver", wait_min, wait_max, wait_avg, 100*wait_avg/delta, wait_std);
        } 
        if (params.verbose) {
            double outgoing_fraction = 1. - ((double) params.recv_per_node) / params.n_recv;
            double volume = sizeof(u64) * N / params.n_send * outgoing_fraction;  // outgoing bytes per node
            char frate[8], nrate[8];
            double delta = wtime() - phase_start;
            human_format(N / params.n_send / delta, frate);
            human_format(volume / delta, nrate);
            printf("phase %d: %.1fs.  %s f/s per process, %sB/s outgoing per node\n", phase, delta, frate, nrate);
        }
    } // phase

    // deal with the results
    std::vector<int> recvcounts(params.size);
    recvcounts[params.rank] = 2 * result.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, params.world_comm);

    std::vector<int> displs(params.size);
    int acc = 0;
    for (int i = 0; i < params.size; i++) {
        displs[i] = acc;
        acc += recvcounts[i];
    }

    std::vector<u64> tmp(acc);
    for (size_t i = 0; i < result.size(); i++) {
        int offset = displs[params.rank] + 2*i;
        auto [x0, x1] = result[i];
        tmp[offset] = x0;
        tmp[offset + 1] = x1;
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT, tmp.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, params.world_comm);
    result.clear();
    for (int i = 0; i < acc; i += 2)
        result.push_back(std::pair(tmp[i], tmp[i+1]));
    if (params.verbose)
        printf("Total: %.1fs\n", wtime() - start);
    return result;
}

}

#endif
