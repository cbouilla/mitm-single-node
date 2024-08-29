#ifndef MITM_NAIVE_MPI_ISEND
#define MITM_NAIVE_MPI_ISEND

#include <vector>
#include <cassert>
#include <cmath>

#include "problem.hpp"
#include "mpi/common.hpp"
#include "dict.hpp"

#include <mpi.h>

/*
 * naive MITM w/ distributed dictionnary.  round-based Alltoallv version.
 */

namespace mitm {

template <bool EXPENSIVE_F, class Problem>
vector<pair<u64, u64>> naive_mpi_claw_search_isend(const Problem &pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    double start = wtime();
    u64 N = 1ull << pb.n;
    vector<pair<u64, u64>> result;
    CompactDict dict((params.role == RECEIVER) ? (1.5 * N) / params.n_recv : 0);

    if (params.verbose) {
        printf("Claw-finding: {0,1}^%d --> {0,1}^%d\n", pb.n, pb.m);
        char hbsize[8], hdsize[8];
        u64 bsize_node = 4 * sizeof(u64) * params.buffer_capacity * params.n_send * params.n_recv / params.n_nodes;
        human_format(bsize_node, hbsize);
        u64 dsize_node = (1.25 * N) / params.n_recv * (sizeof(u64) + sizeof(u32)) * params.recv_per_node;
        human_format(dsize_node, hdsize);
        printf("RAM per node == %sB buffer + %sB dict\n", hbsize, hdsize);
    }

    u64 ncoll = 0;
    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        if (params.verbose)
            printf("Starting phase %d\n", phase);

        double phase_start = wtime();        
        double wait;

        if (params.role == SENDER) {
            SendBuffers sendbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            u64 lo = params.local_rank * N / params.n_send;
            u64 hi = (params.local_rank + 1) * N / params.n_send;
            for (u64 x = lo; x < hi; x++) {
                u64 z = (phase == 0) ? pb.f(x) : pb.g(x);
                u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
                int target = ((int) hash) % params.n_recv;
                if (EXPENSIVE_F)
                    sendbuf.push2(x, z, target);
                else
                    sendbuf.push(x, target);
            }
            sendbuf.flush();

            /* aggregate stats over all senders */
            wait = sendbuf.waiting_time;
        }

        if (params.role == RECEIVER) {
            RecvBuffers recvbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            u64 keys[3 * pb.n];
            while (not recvbuf.complete()) {
                auto ready_buffers = recvbuf.wait();
                for (auto it = ready_buffers.begin(); it != ready_buffers.end(); it++) {
                    auto * buffer = *it;
                    // printf("got buffer! phase=%d, size=%zd\n", phase, buffer->size());
                    for (auto jt = buffer->begin(); jt != buffer->end(); jt++) {
                        u64 x = *jt, z;
                        if (EXPENSIVE_F) {
                            jt++;
                            z = *jt;
                        } else {
                            z = (phase == 0) ? pb.f(x) : pb.g(x);
                        }

                        switch (phase) {
                        case 0:
                            dict.insert(z, x);
                            break;
                        case 1: {
                            // probe dict
                            int nkeys = dict.probe(z, keys);
                            for (int k = 0; k < nkeys; k++) {
                                u64 y = keys[k];
                                if (z != pb.f(y))
                                    continue;    // false positive from truncation in the hash table
                                ncoll += 1;
                                if (pb.is_good_pair(y, x))
                                    result.push_back(pair(y, x));
                            }
                        }
                        }
                    }
                }
            }
            wait = recvbuf.waiting_time;
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
        MPI_Allreduce(MPI_IN_PLACE, &ncoll, 1, MPI_UINT64_T, MPI_SUM, params.world_comm);
        wait_std = std::sqrt(wait_std);
        if (params.local_rank == 0) {
            printf("phase %d %s, wait min %.2fs max %.2fs avg %.2fs (%.1f%%) std %.2fs.\n",
                phase, (params.role == SENDER) ? "sender" : "receiver", wait_min, wait_max, wait_avg, 100*wait_avg/delta, wait_std);
        }
        if (params.verbose) {
            double outgoing_fraction = 1. - ((double) params.recv_per_node) / params.n_recv;
            double volume = sizeof(u64) * N / params.n_nodes * outgoing_fraction;  // outgoing bytes per node
            char frate[8], nrate[8];
            double delta = wtime() - phase_start;
            human_format(N / params.n_send / delta, frate);
            human_format(volume / delta, nrate);
            printf("phase %d: %.1fs.  %s f/s per process, %sB/s outgoing per node. 2^%.2f collisions\n", 
                phase, delta, frate, nrate, std::log2(ncoll));
        }
    } // phase

    if (params.verbose)
        printf("Total: %.1fs\n", wtime() - start);
    
    BCast_result(params, result);
    return result;
}

}

#endif
