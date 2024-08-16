#ifndef MITM
#define MITM

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

class CompactDict {
public:
    const u64 n_slots;     /* How many slots a dictionary have */
    struct __attribute__ ((packed)) entry { u32 k; u64 v; };

    std::vector<struct entry> A;

    CompactDict(u64 n_slots) : n_slots(n_slots)
    {
        A.resize(n_slots, {0xffffffff, 0});
    }

    void insert(u64 key, u64 value)
    {
        u64 h = (key ^ (key >> 32)) % n_slots;
        for (;;) {
            if (A[h].k == 0xffffffff)
                break;
            h += 1;
            if (h == n_slots)
                h = 0;
        }
        A[h].k = key % 0xfffffffb;
        A[h].v = value;
    }

    // return possible values matching this key
    // TODO: replace this by a custom iterator
    int probe(u64 bigkey, u64 keys[])
    {
        u32 key = bigkey % 0xfffffffb;
        u64 h = (bigkey ^ (bigkey >> 32)) % n_slots;
        for (;;) {
            if (A[h].k == 0xffffffff)
                return 0;   // empty-slot, fast path
            if (A[h].k == key)
                break;
            h += 1;
            if (h == n_slots)
                h = 0;
        }
        int nkeys = 0;
        // here, first matching key.
        while (A[h].k == key) {
            keys[nkeys] = A[h].v;
            nkeys += 1;
            h += 1;
            if (h == n_slots)
                h = 0;
        }
        return nkeys;
    }
};


template <class AbstractProblem>
std::vector<std::pair<u64, u64>> naive_mpi_claw_search(AbstractProblem &Pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    double start = wtime();
    MpiCounters ctr;
    u64 N = 1ull << Pb.n;
    std::vector<std::pair<u64, u64>> result;
    CompactDict dict((params.role == RECEIVER) ? (1.5 * N) / params.n_recv : 0);

    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        if (params.verbose)
            printf("Starting phase %d\n", phase);
        double phase_start = wtime();
        ctr.reset();
        
        if (params.role == SENDER) {
            BaseSendBuffers sendbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            u64 lo = params.local_rank * N / params.n_send;
            u64 hi = (params.local_rank + 1) * N / params.n_send;
            for (u64 x = lo; x < hi; x++) {
                u64 z = (phase == 0) ? Pb.f(x) : Pb.g(x);
                u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
                int target = ((int) hash) % params.n_recv;
                sendbuf.push(x, target, ctr);
            }
            sendbuf.flush(ctr);

            char frate[8], nrate[8];
            double delta = wtime() - phase_start;
            human_format(N / params.n_send / delta, frate);
            human_format(ctr.bytes_sent / delta, nrate);
            printf("phase %d, sender %d, wait %.3fs (%.1f%%), %s f/s, %sB/s\n", 
                phase, params.local_rank, ctr.send_wait, 100*ctr.send_wait/delta, frate, nrate);
        }

        if (params.role == RECEIVER) {
            BaseRecvBuffers recvbuf(params.inter_comm, TAG_POINTS, params.buffer_capacity);
            while (not recvbuf.complete()) {
                auto ready_buffers = recvbuf.wait(ctr);
                for (auto it = ready_buffers.begin(); it != ready_buffers.end(); it++) {
                    auto * buffer = *it;
                    // printf("got buffer! phase=%d, size=%zd\n", phase, buffer->size());
                    for (auto jt = buffer->begin(); jt != buffer->end(); jt++) {
                        u64 x = *jt;
                        switch (phase) {
                        case 0: 
                            dict.insert(Pb.f(x), x);
                            break;
                        case 1: {
                            // probe dict
                            u64 z = Pb.g(x);
                            u64 keys[3* Pb.n];
                            int nkeys = dict.probe(z, keys);
                            for (int k = 0; k < nkeys; k++) {
                                u64 y = keys[k];
                                if (z != Pb.f(y)) {
                                    ctr.collision_failure();
                                    continue;    // false positive from truncation in the hash table
                                }
                                ctr.found_collision();
                                if (Pb.is_good_pair(y, x)) {
                                    printf("\nfound golden collision !!!\n");
                                    result.push_back(std::pair(y, x));
                                }
                            }
                        }
                        }
                    }
                }
            }
            char frate[8];
            double delta = wtime() - phase_start;
            human_format(N / params.n_recv / delta, frate);
            printf("phase %d, sender %d, wait %.3fs (%.1f%%), %s f/s\n", 
                phase, params.local_rank, ctr.recv_wait, 100*ctr.recv_wait/delta, frate);
            // deal with result
        } // RECEIVER
        if (params.verbose)
            printf("Phase: %.1fs\n", wtime() - phase_start);
    } // phase
    if (params.verbose)
        printf("Total: %.1fs\n", wtime() - start);
    return result;
}

}

#endif
