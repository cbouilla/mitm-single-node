#ifndef MITM
#define MITM

#include <optional>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cmath>

#include "common.hpp"
#include "AbstractCollisionProblem.hpp"

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
    std::pair<u64, u64> probe(u64 bigkey)
    {
        u32 key = bigkey % 0xfffffffb;
        u64 h = (bigkey ^ (bigkey >> 32)) % n_slots;
        for (;;) {
            if (A[h].k == 0xffffffff)
                return std::pair(0, 0);   // empty-slot, fast path
            if (A[h].k == key)
                break;
            h += 1;
            if (h == n_slots)
                h = 0;
        }
        // here, first matching key.
        u64 lo = h;
        while (A[h].k == key) {
            h += 1;
            if (h == n_slots)
                h = 0;            
        }
        return std::pair(lo, h);
    }
};

template <class AbstractProblem>
std::vector<std::pair<u64, u64>> naive_mpi_claw_search(AbstractProblem &Pb)
{
    static_assert(std::is_base_of<AbstractClawProblem, AbstractProblem>::value,
        "problem not derived from mitm::AbstractClawProblem");
  
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start = wtime();
    u64 N = 1ull << Pb.n;
    CompactDict dict(2 * N/size);
    std::vector<std::pair<u64, u64>> result;

    const u64 alpha = 10000;    // expected #values received in each round by each process
    const u64 K = alpha*size;   // total values generated in each round
    // entries for each target follows binomial law (#values, 1/size);

    // generate K values :
    // each bucket <= alpha + sqrt(210 * alpha)   with proba 2^{-100}
    int limit = alpha + std::sqrt(210 * alpha);
    if (rank == 0)
        printf("limit=%d\n", limit);

    u64 probe_false_pos = 0;
    std::vector<u64> sendbuffer(size * limit);
    std::vector<u64> recvbuffer(size * limit);
    std::vector<int> sendcounts(size);
    std::vector<int> recvcounts(size);
    std::vector<int> displs(size);
    for (int i = 0; i < size; i++)
        displs[i] = limit * i;

    if (rank == 0) {
        char hsize[8];
        human_format(16 * size * limit, hsize);
        printf("Buffer size (per process) = %s\n", hsize);
    }
    double wait = 0;

    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        if (rank == 0)
            printf("Starting phase %d\n", phase);
        double phase_start = wtime();

        const u64 nrounds = (N + K - 1) / K;
        for (u64 round = 0; round < nrounds; round ++) {
            // reset sendcounts
            for (int i = 0; i < size; i++)
                sendcounts[i] = 0;

            // round = [N * round / nrounds : N * (round + 1) / nrounds]
            u64 round_lo = N * round / nrounds;
            u64 round_hi = N * (round + 1) / nrounds;
            u64 round_size = round_hi - round_lo;
            u64 process_lo = round_lo + rank * round_size / size;
            u64 process_hi = round_lo + (rank + 1) * round_size / size;
            for (u64 x = process_lo; x < process_hi; x++) {
                u64 z = (phase == 0) ? Pb.f(x) : Pb.g(x);
                u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
                int target = ((int) hash) % size;
                assert(sendcounts[target] < limit);
                u64 offset = limit * target + sendcounts[target];
                sendbuffer[offset] = x;
                sendcounts[target] += 1;
            }
    
        double start_comm = wtime();

        // exchange buffer sizes;
        MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        // exchange data / TODO: IAlltoallv
        MPI_Alltoallv(sendbuffer.data(), sendcounts.data(), displs.data(), MPI_UINT64_T, 
                      recvbuffer.data(), recvcounts.data(), displs.data(), MPI_UINT64_T, MPI_COMM_WORLD);

        wait += wtime() - start_comm;

        for (int i = 0; i < size; i++)
            for (int j = 0; j < recvcounts[i]; j++) {
                u64 x = recvbuffer[i * limit + j];
                if (phase == 0) {
                    // inject data into dict
                    u64 z = Pb.f(x);
                    dict.insert(z, x);
                } else {
                    // probe dict
                    u64 z = Pb.g(x);
                    auto [lo, hi] = dict.probe(z);
                    for (u64  k = lo; k < hi; k++) {
                        u64 y = dict.A[k].v;
                        if (z != Pb.f(y)) {
                            probe_false_pos += 1;
                            continue;    // false positive from truncation in the hash table
                        }
                        if (Pb.is_good_pair(y, x)) {
                            printf("\nfound golden collision !!!\n");
                            result.push_back(std::pair(y, x));
                        }
                    }
                }
            }

        if (rank == 0) {
            char frate[8], nrate[8];
            double delta = wtime() - start;
            human_format(K * (round + 1) / delta, frate);
            human_format(8 * K * (round + 1) / delta, nrate);
            printf("\rRound %" PRId64 " / %" PRId64 ".  Wait/round = %.3fs.  %s f()/s.  Net=%sB/s", 
                    round, nrounds, wait/(1+round), frate, nrate);

            fflush(stdout);
        }
    }
    if (rank == 0)
        printf("\nPhase: %.1fs\n", wtime() - phase_start);
    }
    return result;
}

}

#endif
