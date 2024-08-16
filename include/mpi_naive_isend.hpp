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

/* try to iterate for 1s. Return #it/s */
template<typename AbstractProblem>
static void naive_benchmark(const AbstractProblem& Pb, MpiParameters &params)
{
    MPI_Barrier(params.world_comm);

    u64 N = 1ull << 30; 
    double start = wtime();
    u64 count = 0;
    for (u64 x = 0; x < N; x++) {
        u64 z = (x & 1) ? Pb.f(x) : Pb.g(x);
        u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
        int target = ((int) hash) % params.n_recv;
        if (target == 0)
            count += 1;
    }
    double rate = wtime() - start;
    double rate_min = rate;
    double rate_max = rate;
    double rate_avg = rate;
    MPI_Allreduce(MPI_IN_PLACE, &rate_min, 1, MPI_DOUBLE, MPI_MIN, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_max, 1, MPI_DOUBLE, MPI_MAX, params.world_comm);
    MPI_Allreduce(MPI_IN_PLACE, &rate_avg, 1, MPI_DOUBLE, MPI_SUM, params.world_comm);
    rate_avg /= params.size;
    double rate_std = (rate - rate_avg) * (rate - rate_avg);
    MPI_Allreduce(MPI_IN_PLACE, &rate_std, 1, MPI_DOUBLE, MPI_SUM, params.local_comm);
    rate_std = std::sqrt(rate_std);
    if (params.rank == 0) {
        char hmin[8], hmax[8], havg[8], hstd[8];
        human_format(rate_min, hmin);
        human_format(rate_max, hmax);
        human_format(rate_avg, havg);
        human_format(rate_std, hstd);
        printf("Benchmark. f/s: min %s max %s avg %s std %s (aggregated over %d processes)\n", hmin, hmax, havg, hstd, params.size);
    }
}


template <bool EXPENSIVE_F, class AbstractProblem>
std::vector<std::pair<u64, u64>> naive_mpi_claw_search_isend(const AbstractProblem &Pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    double start = wtime();
    MpiCounters ctr;
    u64 N = 1ull << Pb.n;
    std::vector<std::pair<u64, u64>> result;
    CompactDict dict((params.role == RECEIVER) ? (1.5 * N) / params.n_recv : 0);

    if (params.verbose) {
        char hbsize[8], hdsize[8];
        u64 bsize_node = 4 * sizeof(u64) * params.buffer_capacity * params.n_send * params.n_recv / params.n_nodes;
        human_format(bsize_node, hbsize);
        u64 dsize_node = (1.25 * N) / params.n_recv * (sizeof(u64) + sizeof(u32)) * params.recv_per_node;
        human_format(dsize_node, hdsize);
        printf("RAM per node == %sB buffer + %sB dict\n", hbsize, hdsize);
    }

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
                if (EXPENSIVE_F)
                    sendbuf.push2(x, z, target, ctr);
                else
                    sendbuf.push(x, target, ctr);
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
                        u64 x = *jt, z;
                        if (EXPENSIVE_F) {
                            jt++;
                            z = *jt;
                        } else {
                            z = (phase == 0) ? Pb.f(x) : Pb.g(x);
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
            double volume = sizeof(u64) * N / params.n_nodes * outgoing_fraction;  // outgoing bytes per node
            char frate[8], nrate[8];
            double delta = wtime() - phase_start;
            human_format(N / params.n_send / delta, frate);
            human_format(volume / delta, nrate);
            printf("phase %d: %.1fs.  %s f/s per process, %sB/s outgoing per node\n", phase, delta, frate, nrate);
        }
    } // phase

    if (params.verbose)
        printf("Total: %.1fs\n", wtime() - start);
    
    BCast_result(params, result);
    return result;
}

}

#endif
