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


/* Manages send buffers for a collection of receiver processes, with double-buffering */
class SendBuffers {
    using Buffer = std::vector<u64>;
public:
    double wait = 0.;
    u64 bytes_sent = 0;

private:
    const size_t capacity;
    int n;

    std::vector<Buffer> ready;
    std::vector<Buffer> outgoing;
    std::vector<MPI_Request> request;   /* for the OUTGOING buffers */

    /* wait until the i-th passive buffer has been fully sent */
    void wait_send(int i)
    {
        double start = wtime();
        MPI_Wait(&request[i], MPI_STATUS_IGNORE);
        wait += wtime() - start;
        outgoing[i].clear();
    }

    /* initiate transmission of the i-th OUTGOING buffer */
    void start_send(int i)
    {
        if (outgoing[i].size() == 0)  // do NOT send empty buffers. These are interpreted as "I am done"
            return;
        MPI_Issend(outgoing[i].data(), outgoing[i].size(), MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &request[i]);
        bytes_sent += outgoing[i].size() * sizeof(u64);
    }

public:
    SendBuffers(size_t capacity) : capacity(capacity)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &n);
        ready.resize(n);
        outgoing.resize(n);
        request.resize(n, MPI_REQUEST_NULL);
        for (int i = 0; i < n; i++) {
            ready[i].reserve(capacity);
            outgoing[i].reserve(capacity);
        }
    }

    ~SendBuffers()
    {
        for (int i = 0; i < n; i++)
            assert(request[i] == MPI_REQUEST_NULL);
    }

    /* add a new item to the send buffer. Send if necessary */
    void push(u64 x, int rank)
    {
        if (ready[rank].size() == capacity) {
            /* ready buffer is full. */
            wait_send(rank);      // finish sending the outgoing buffer
            std::swap(ready[rank], outgoing[rank]);
            start_send(rank);     // start sending
        }
        ready[rank].push_back(x);
    }

    /* send and empty all buffers, even if they are incomplete */
    void flush()
    {
        // finish sending all the outgoing buffers
        double start = wtime();
        MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

        // send all the (incomplete) ready buffers
        for (int i = 0; i < n; i++) {
            std::swap(ready[i], outgoing[i]);
            start_send(i);
        }
        MPI_Waitall(n, request.data(), MPI_STATUSES_IGNORE);

        // finally tell all receivers that we are done
        for (int i = 0; i < n; i++)
            MPI_Send(NULL, 0, MPI_UINT64_T, i, 0, MPI_COMM_WORLD);
        wait += wtime() - start;
    }
};



/* Manage reception buffers for a collection of sender processes, with double-buffering */
class RecvBuffers {
    using Buffer = std::vector<u64>;
public:
    double wait = 0.;
    u64 bytes_sent = 0;

private:
    const size_t capacity;
    int n;

    std::vector<Buffer> ready;                 // buffers containing points ready to be processed 
    std::vector<Buffer> incoming;              // buffers waiting for incoming data
    std::vector<MPI_Request> request;
    int n_active_senders;                      // # active senders

    /* initiate reception for a specific sender */
    void listen_sender(int i)
    {
        incoming[i].resize(capacity);
        MPI_Irecv(incoming[i].data(), capacity, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &request[i]);
    }

public:
    RecvBuffers(size_t capacity) : capacity(capacity)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &n);
        ready.resize(n);
        incoming.resize(n);
        request.resize(n, MPI_REQUEST_NULL);
        for (int i = 0; i < n; i++) {
            ready[i].reserve(capacity);
            incoming[i].reserve(capacity);
            listen_sender(i);
        }
        n_active_senders = n;

    }

    ~RecvBuffers()
    {
        assert(n_active_senders == 0);
        for (int i = 0; i < n; i++)
            assert(request[i] == MPI_REQUEST_NULL);
    }

    /* return true when all senders are done */
    bool complete()
    {
        return (n_active_senders == 0);
    }

    /* 
     * return all the ready buffers (since last call).
     * This may destroy the content of all previously ready buffers, so that they have to be processed first
     */
    std::vector<Buffer *> probe()
    {
        std::vector<Buffer *> result;
        int n_done;
        std::vector<int> rank_done(n);
        std::vector<MPI_Status> statuses(n);
        MPI_Testsome(n, request.data(), &n_done, rank_done.data(), statuses.data());
        
        for (int i = 0; i < n_done; i++) {
            int j = rank_done[i];
            std::swap(incoming[j], ready[j]);
            int count;
            MPI_Get_count(&statuses[i], MPI_UINT64_T, &count);
            if (count == 0) {
                n_active_senders -= 1;
            } else {
                ready[j].resize(count);         // matching message size
                result.push_back(&ready[j]);
                listen_sender(j);
            }
        }
        return result;
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
    CompactDict dict((1.25 * N) / size);
    std::vector<std::pair<u64, u64>> result;

    const u64 alpha = 2048;    // expected #values received in each round by each process
    const u64 K = alpha*size;   // total values generated in each round
    // entries for each target follows binomial law (#values, 1/size);

    // generate K values :
    // each bucket <= alpha + sqrt(210 * alpha)   with proba 2^{-100}
    int limit = alpha + std::sqrt(120 * alpha);
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

    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        if (rank == 0)
            printf("Starting phase %d\n", phase);
        double phase_start = wtime();
	double last_display = phase_start;
	double wait = 0;

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

	double now = wtime();
        if (rank == 0 && now - last_display >= 0.5) {
            char frate[8], nrate[8];
            double delta = now - phase_start;
	    last_display = now;
            human_format(K * (round + 1) / delta, frate);
            human_format(8 * K * (round + 1) / delta, nrate);
            printf("Round %" PRId64 " / %" PRId64 ".  Wait/round = %.3fs (%.1f%%).  %s f()/s.  Net=%sB/s\n",
		   round, nrounds, wait/(1+round), 100.*wait / delta, frate, nrate);

            fflush(stdout);
        }
    }
    if (rank == 0)
        printf("Phase: %.1fs\n", wtime() - phase_start);
    }
    if (rank == 0)
        printf("Total: %.1fs\n", wtime() - start);
    return result;
}

}

#endif
