#ifndef MITM_MPI_DIRECT
#define MITM_MPI_DIRECT

#include <vector>
#include <queue>
#include <cassert>
#include <cmath>

#include <mpi.h>
#include <omp.h>

#include "problem.hpp"
#include "mpi/common.hpp"
#include "dict.hpp"

using std::atomic;

/*
 * direct MITM with distributed dictionnary. MPI + OpenMP version
 */

namespace mitm {

enum tags {TAG_POINTS, TAG_MANAGEMENT};
enum service {DONE};

class LockfreeStack;

class Buffer : public vector<u64> {
public:
    Buffer* stack_next;
    LockfreeStack *stack;
    bool busy = false;
};

template <class T>
bool CAS(std::atomic<T> &target, T expected, T desired)
{
    return target.compare_exchange_strong(expected, desired);
    // bool ok = false;
    // #pragma omp atomic compare capture
    // { 
    //     ok = target == expected; if (ok) { target = desired; }
    // }
    // return ok;
}

class LockfreeStack {
private:
    atomic<Buffer *> head;
    atomic<size_t> n;
public:
    LockfreeStack()
    {
        n = 0;
        head = nullptr;
    };

    void push(Buffer *bufptr)
    {
        assert(bufptr != nullptr);
        assert(bufptr->stack == nullptr);
        bufptr->stack = this;
        for (;;) {
            bufptr->stack_next = head.load();
            if (CAS(head, bufptr->stack_next, bufptr)) {
                n += 1;
                return;
            }
        }
    }

    Buffer * try_pop()
    {
        Buffer *bufptr;
        for (;;) {
            bufptr = head.load();
            if (bufptr == nullptr)
                return nullptr;
            if (CAS(head, bufptr, bufptr->stack_next))
                break;
        }
        assert(bufptr->stack == this);
        bufptr->stack = nullptr;
        n -= 1;
        return bufptr;
    }

    bool empty()
    {
        return (head.load() == nullptr);
    }

    size_t size()
    {
        return n.load();
    }
};

/* Manage a collection of buffers.  The n vectors are **reserved** (but not **resized**) to capacity */
class FreeList {
private:
    LockfreeStack free;
    vector<Buffer> all;

public:
    void setup(size_t capacity, size_t n)
    {
        all.resize(n);
        for (size_t i = 0; i < n; i++) {
            all[i].reserve(capacity);
            assert(not all[i].busy);
            free.push(&all[i]);
        }
    }

    void release(Buffer *bufptr)
    {
        assert(bufptr->busy);
        bufptr->busy = 0;
        free.push(bufptr);
    }

    Buffer * try_acquire()
    {
        Buffer *bufptr = free.try_pop();
        if (bufptr == nullptr)
            return nullptr;
        assert(not bufptr->busy);
        bufptr->busy = 1;
        return bufptr;
    }

    size_t size()
    {
        return free.size();
    }

    size_t capacity()
    {
        return all.size();
    }

    bool full()
    {
        return (size() == capacity());
    }

    void clear()
    {
        assert(full());
    }
};

enum sendstate {READY, SIZE, SEND};
enum recvstate {PASSIVE, RECV, DISCHARGE, GRACE};

template <class Problem>
class ClawSearchProcess {
private:

    struct recv_buffer {
        vector<u64> payload;
        atomic<int> state;
        atomic<size_t> lo;
        atomic<size_t> hi;
        atomic<size_t> done;
    };

    FreeList freelist;                             // buffers exchanged between threads
    atomic<int> n_workers_started = 0;                     
    atomic<int> n_workers_done = 0;
    vector<u64> bytes_sent;
    double start;

    // communication between workers and coordinator
    LockfreeStack *outgoing;                       // buffers waiting to be sent, one queue per MPI rank
    atomic<bool> coordinator_done;

    MpiParameters &params;
    MPI_Comm comm;                                  // shorthand for params.comm
    u64 chunksize;                                  // shorthand for params.chunksize
    size_t bufratio;                                // shorthand for params.inter_buffer_capacity / params.inter_buffer_capacity;
    int sendstate;
    struct recv_buffer recvbuf[2];
    MPI_Request req;

    int phase = 0;

    // producer state
    u64 N;
    u64 sender_lo, sender_hi;
    
    /*
     * coordination between threads
     * ============================
     * - A worker thread is DONE when it will no longer enqueue an outgoing buffer.
     * - When all the workers in a host are DONE, and all the outgoing buffers have been completely sent, then this host is DONE.
     *
     */ 

    CompactDict dict;
    
    /*
    void process_service_message(int source, Buffer *bufptr)
    {
        Buffer &buf = *bufptr;
        assert(buf.size() > 0);
        switch(buf[0]) {
        case DONE:
            deactivate_host();
            break;
        }
    }*/

    void worker_start()
    {
        n_workers_started += 1;
    }

    void worker_done()
    {
        n_workers_done += 1;
        // println("mpi rank {}: {} / {} Worker done!", params.mpi_rank, n_workers_done.load(), n_workers_started.load());
    }

    bool all_workers_done()
    {
        return (n_workers_started.load() > 0) && (n_workers_started.load() == n_workers_done.load());
    }

    bool finished()
    {
        return coordinator_done.load() or (params.mpi_size == 1);
    }

    void verbosity()
    {
        // worker progress ; "hash rate" ; network (MB/s, busy send, busy recv); mem (available buffers)
        u64 range_size = N / params.mpi_size;
        u64 remaining = sender_hi - sender_lo;
        u64 done = range_size - remaining;
        double progress = 100. * done / range_size;
        char hfrate[8], hnrate[8];
        double delta = wtime() - start;
        human_format(done / delta, hfrate);
        u64 total_bsent = 0;
        for (int i = 0; i < params.mpi_size; i++)
            total_bsent += bytes_sent[i];
        human_format(total_bsent / delta, hnrate);
        int in_flight = freelist.capacity() - freelist.size() - (params.n_threads - 1 - n_workers_done) * params.mpi_size;
        println("\rDone: {:.1f}%. {} f/s. net: {}B/s. {} in-flight intra buffers. Sendsate={}, recvstate=({}; {})",
            progress, hfrate, hnrate, in_flight, sendstate, recvbuf[0].state.load(), recvbuf[1].state.load());
        std::fflush(stdout);
    }

    int get_recvbuf(int desired)
    {
        for (int i = 0; i < 2; i++)
            if (recvbuf[i].state.load() == desired)
                return i;
        return -1;
    }

public:
    const Problem &pb;
    vector<pair<u64, u64>> result;
    atomic<u64> ncoll = 0;


    ClawSearchProcess(const Problem &pb, MpiParameters &_params) : 
        params(_params), comm(params.comm), chunksize(params.chunksize), pb(pb)
    {
        assert(params.n_threads != 0);
        bufratio = params.inter_buffer_capacity / params.intra_buffer_capacity;
        u64 intra_inflight = params.mpi_size * params.n_threads + params.n_recv_buffers * bufratio;
        u64 n_intra_buffers = params.send_buffers_ratio * intra_inflight + params.n_slack_buffers;
        N = 1ull << pb.n;
        ncoll = 0;

        freelist.setup(params.intra_buffer_capacity, n_intra_buffers);
        dict.set_size((params.dict_capacity_ratio * N) / params.mpi_size);
        outgoing = new LockfreeStack[params.mpi_size];
        u64 limit = params.inter_buffer_capacity;
        for (int i = 0; i < 2; i++)
            recvbuf[i].payload.resize(limit * params.mpi_size);

        if (params.verbose) {
            char hbsize[8], hdsize[8];
            u64 intrasize = n_intra_buffers * params.intra_buffer_capacity * sizeof(u64);
            u64 intersize = 2 * params.mpi_size * params.inter_buffer_capacity;
            human_format(intrasize + intersize, hbsize);
            human_format(dict.nbytes(), hdsize);
            println("RAM per node == {}B buffer () +{}B dict ({} slots)", hbsize, hdsize, dict.n_slots);
        }
    }

    ~ClawSearchProcess()
    {
        delete outgoing;
    }

    void start_phase(int i)
    {
        // application-specific
        sender_lo = params.mpi_rank * N / params.mpi_size;
        sender_hi = (params.mpi_rank + 1) * N / params.mpi_size;
        std::string name[2] = {"fill", "probe"}; 
        if (params.verbose)
            println("starting phase {} ({})", i, name[i]);
        
        // infrastructure
        freelist.clear();
        phase = i;
        n_workers_started = 0;
        n_workers_done = 0;
        start = wtime();
        bytes_sent.clear();
        bytes_sent.resize((size_t) params.mpi_size);
        for (int i = 0; i < params.mpi_size; i++)
            assert(outgoing[i].empty());
    }

    void coordinator()
    {
        u64 limit = params.inter_buffer_capacity;
        vector<u64> sendbuf(limit * params.mpi_size);
        
        // ??????? prepare a buffer for the MPI_Bsends;
        u64 mpi_buffer_size = params.mpi_size * (8 + MPI_BSEND_OVERHEAD);
        u8 mpi_buffer[mpi_buffer_size];
        MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);

        int backoff = params.min_backoff;
        double last_verbose = start;
        req = MPI_REQUEST_NULL;
        sendstate = READY;
        recvbuf[0].state = PASSIVE;
        recvbuf[1].state = PASSIVE;
        vector<int> sendcounts(params.mpi_size);
        vector<int> recvcounts(params.mpi_size);
        vector<int> senddispls(params.mpi_size);
        vector<int> recvdispls(params.mpi_size);
        for (int i = 0; i < params.mpi_size; i++)
            senddispls[i] = i * limit;
        int receiver = -1;                               // track the receiving buffer

        coordinator_done = 0;
        bool signaled_workers_done = 0;
        while (not coordinator_done.load()) {                   // loop until the action has died down
            bool action = 0;

            bool flushing = all_workers_done();
            if (flushing and not signaled_workers_done) {
                signaled_workers_done = 1;
                println("MPI rank {}, all workers done ({:.1f}s)", params.mpi_rank, wtime() - start);
            }

            // locate a passive recv buffer
            int passive_recvbuf = get_recvbuf(PASSIVE);
            
            // decide if a new round must start
            bool trigger = 1;
            for (int i = 0; i < params.mpi_size; i++)
                if (i != params.mpi_rank and (outgoing[i].size() < bufratio) and not flushing)
                    trigger = 0;
            if (sendstate != READY or passive_recvbuf < 0)
                trigger = 0;

            if (trigger) {
                assert(req == MPI_REQUEST_NULL);
                action = 1;
                bool done = flushing and freelist.full();
                // copy blocks from the pending queues to the send buffer
                for (int i = 0; i < params.mpi_size; i++) {
                    sendcounts[i] = (done ? -1 : 0);
                    while (Buffer *bufptr = outgoing[i].try_pop()) {
                        Buffer &buf = *bufptr;
                        if (sendcounts[i] + buf.size() > limit) {
                            outgoing[i].push(bufptr);  // too large, put it back
                            break;
                        }
                        // copy the buffer and release it
                        u64 start = senddispls[i] + sendcounts[i];
                        for (u64 k = 0; k < buf.size(); k++)
                            sendbuf[start + k] = buf[k];
                        sendcounts[i] += buf.size();
                        freelist.release(bufptr);
                    }
                    assert(sendcounts[i] <= (int) limit);
                }
                sendstate = SIZE;
                MPI_Ialltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD, &req);
                // println("Send to {} started ({} items) [flushing={} / active={}]", i, sendbuf[i].size(), flushing, n_active_sends);
            }
    
            if (sendstate != READY) {
                int flag;
                MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
                if (flag) {
                    if (sendstate == SEND) {
                        // make the data available to the workers
                        assert(receiver >= 0);
                        assert(recvbuf[receiver].state.load() == RECV);
                        recvbuf[receiver].state = DISCHARGE;
                        receiver = -1;
                        sendstate = READY;
                        for (int i = 0; i < params.mpi_size; i++)
                            bytes_sent[i] += sizeof(u64) * sendcounts[i];

                    } else if (sendstate == SIZE) {
                        assert(passive_recvbuf >= 0);
                        assert(recvbuf[passive_recvbuf].state.load() == PASSIVE);
                        // detect potential termination
                        bool finished_sending = 1;
                        for (int i = 0; i < params.mpi_size; i++) {
                            if (sendcounts[i] < 0)
                                sendcounts[i] = 0;
                            if (recvcounts[i] < 0)
                                recvcounts[i] = 0;
                            else
                                finished_sending = 0;
                        }

                        if (finished_sending and recvbuf[0].state.load() == PASSIVE and recvbuf[1].state.load() == PASSIVE) {
                            coordinator_done = 1;
                            continue;                              // exit the coordinator main loop
                        }
                        // prepare recvdispl to stack all incoming data
                        int acc = 0;
                        for (int i = 0; i < params.mpi_size; i++) {
                            recvdispls[i] = acc;
                            acc += recvcounts[i];
                        }
                        recvbuf[passive_recvbuf].done = 0;
                        recvbuf[passive_recvbuf].lo = 0;
                        recvbuf[passive_recvbuf].hi = acc;
                        // maybe Ialltoall if the sizes match?
                        MPI_Ialltoallv(sendbuf.data(), sendcounts.data(), senddispls.data(), MPI_UINT64_T, 
                            recvbuf[passive_recvbuf].payload.data(), recvcounts.data(), recvdispls.data(), MPI_UINT64_T, 
                            MPI_COMM_WORLD, &req);
                        sendstate = SEND;
                        receiver = passive_recvbuf;
                        recvbuf[receiver].state = RECV;
                    }
                }
            }

            if (params.verbose and wtime() >= last_verbose + 1) {
                verbosity();
                last_verbose = wtime();
            }

            if (action) {
                backoff = params.min_backoff;
            } else {
                backoff = std::min(params.max_backoff, 2*backoff + 1);
                wait(backoff);
            }
        }
        print("MPI rank {}, coordinator main loop done", params.mpi_rank);

        assert(req == MPI_REQUEST_NULL);
        assert(freelist.full());
        for (int i = 0; i < params.mpi_size; i++)
            assert(outgoing[i].empty());
        assert(recvbuf[0].state.load() == PASSIVE);
        assert(recvbuf[1].state.load() == PASSIVE); 

        for (int i = 0; i < params.mpi_size; i++) {
            MPI_Barrier(params.comm);
            if (i == params.mpi_rank) {
                print("MPI rank {} : ", i);
                for (int j = 0; j < params.mpi_size; j++)
                    if (j != i) {
                        char hbs[8];
                        human_format(bytes_sent[j], hbs);
                        print(" {}B --> {} / ", hbs, j);
                    }
                println("");
            }
        }

        void *foo;
        int bar;
        MPI_Buffer_detach(&foo, &bar);
    }

    /*************** handling of coordinator --> workers buffers *****************/

    tuple<int, size_t, size_t> grab_incoming_chunk()
    {
        tuple<int, size_t, size_t> result;
        #pragma omp critical(recvbuf)
        {
            int i = get_recvbuf(DISCHARGE);
            if (i < 0)
                result = tuple(-1, 0, 0);
            else {
                u64 lo = recvbuf[i].lo;
                u64 next = lo + params.intra_buffer_capacity;
                if (next >= recvbuf[i].hi) {
                    next = recvbuf[i].hi;
                    recvbuf[i].state = GRACE;
                } else {
                    recvbuf[i].lo = next;
                }
                result = tuple(i, lo, next);
            }
        }
        return result;
    }


    void poll_incoming()
    {
        for (;;) {
            auto [i, lo, hi] = grab_incoming_chunk();
            if (i < 0)
                return;
            process_incoming_buffer(recvbuf[i].payload.data(), lo, hi);
            #pragma omp critical(recvbuf)
            {
                recvbuf[i].done += hi - lo;
                if (recvbuf[i].done.load() == recvbuf[i].hi) {      // the atomic<> on done is not necessary
                    assert(recvbuf[i].state == GRACE);
                    recvbuf[i].state = PASSIVE;
                }
            }
        }
    }

    /*************** handling of workers --> coordinator buffers *****************/

    Buffer * get_worker_buffer()
    {
        for (;;) {
            poll_incoming();
            Buffer *bufptr = freelist.try_acquire();
            if (bufptr != nullptr) {
                bufptr->clear();
                return bufptr;
            }
        }
    }

    void send_buffer(Buffer *bufptr, int rank)
    {
        if (rank == params.mpi_rank) {
            // fast-track: the current thread deals with it right now
            process_incoming_buffer(bufptr->data(), 0, bufptr->size());
            freelist.release(bufptr);
        } else {
            // slow track: the coordinator deals with it
            outgoing[rank].push(bufptr);
        }
    }

    /********************* application-specific code *********************************/
    
    // this is specific to the direct mitm
    void process_incoming_buffer(u64 *buf, u64 lo, u64 hi)
    {
        int maxkeys = 3 * pb.n;
        u64 keys[maxkeys];
        
        for (u64 i = lo; i < hi; i += 2) {
            u64 x = buf[i];
            u64 z = buf[i + 1];
            switch (phase) {
            case 0: {
                dict.insert(z, x);
                break;
            }
            case 1: {
                int nkeys = dict.probe(z, maxkeys, keys);
                if (nkeys > maxkeys)
                    println("ERROR: too many keys");
                for (int k = 0; k < nkeys; k++) {
                    u64 y = keys[k];
                    if (z != pb.f(y))
                        continue;    // false positive from truncation in the hash table
                    ncoll += 1;
                    if (pb.is_good_pair(y, x)) {
                        println("!! solution found !!");
                        #pragma omp critical
                        result.push_back(pair(y, x));
                    }
                }
            }
            }
        }
    }

    pair<u64, u64> grab_work()
    {
        u64 local;
        #pragma omp atomic capture 
        { local = sender_lo; sender_lo += chunksize; }
        return pair(local, std::min(local + chunksize, sender_hi));
    }

    // many such thread per node. This is specific to the direct mitm
    void worker()
    {
        vector<Buffer *> buffers;
        for (int i = 0; i < params.mpi_size; i++) {
            Buffer *bufptr = get_worker_buffer();
            buffers.push_back(bufptr);
        }

        worker_start();

        for (;;) {
            poll_incoming();
            auto [chunk_lo, chunk_hi] = grab_work();
            if (chunk_lo >= chunk_hi)
                break;                                // no more work
            for (u64 x = chunk_lo; x < chunk_hi; x++) {
                    u64 z = (phase == 0) ? pb.f(x) : pb.g(x);
                    u64 hash = (z * 0xdeadbeef) % 0x7fffffff;
                    int i = ((int) hash) % params.mpi_size;
                    if (buffers[i]->size() == buffers[i]->capacity()) {
                        send_buffer(buffers[i], i);
                        buffers[i] = get_worker_buffer();
                    }
                    buffers[i]->push_back(x);
                    buffers[i]->push_back(z);
            }
        }
        // flush output buffers
        for (size_t i = 0; i < buffers.size(); i++)
            send_buffer(buffers[i], i);
        
        worker_done();
        
        // process incoming buffers while the whole computation is not finished
        while (not finished())
            poll_incoming();
    }
};


template <class Problem>
vector<pair<u64, u64>> mpi_direct_claw_search(const Problem &pb, MpiParameters &params)
{
    static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
        "problem not derived from mitm::AbstractClawProblem");

    if (params.verbose) {
        println("Claw-finding: {{0,1}}^{} --> {{0,1}}^{}", pb.n, pb.m);
        if (params.mpi_size > 1)
            println("Distributed mode: running on {} MPI processes", params.mpi_size);
        else
            println("Single-host mode (OpenMP only)");
    }

    if (params.n_threads <= 0) {
        params.n_threads = omp_get_max_threads();
        if (params.verbose)
            println("Autodetect: using {} threads", params.n_threads);
    }

    if (params.max_concurrent_send <= 0) 
        params.max_concurrent_send = 2 * params.mpi_size;

    double start = wtime();
    ClawSearchProcess proc(pb, params);    

    for (int phase = 0; phase < 2; phase++) {
        // phase 0 == fill the dict with f()
        // phase 1 == probe the dict with g()
        
        proc.start_phase(phase);
        // double phase_start = wtime();        

        #pragma omp parallel
        {
            assert(params.mpi_size == 1 or omp_get_num_threads() > 1);
            int tid = omp_get_thread_num();
            if (tid == 0 and params.mpi_size > 1)
                proc.coordinator();
            else
                proc.worker();
        }

        #if 0
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
        #endif
    } // phase

    if (params.verbose)
        println("Total: {:.1f}s", wtime() - start);
    
    BCast_result(params, proc.result);
    return proc.result;
};

}

#endif


// 4xgrvingt

// OMP_NUM_THREADS=32 mpirun --map-by ppr:1:node --mca pml ^ucx examples/mpi_double_speck64_direct --n 28
// ---> 332.6s

// mpirun --map-by ppr:1:socket --mca pml ^ucx examples/mpi_double_speck64_direct --n 28   (2x16 threads)
// ---> 490s
