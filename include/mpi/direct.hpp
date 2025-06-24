#ifndef MITM_MPI_DIRECT
#define MITM_MPI_DIRECT

#include <vector>
#include <cassert>
#include <cmath>

#include <mpi.h>
#include <omp.h>

#include "problem.hpp"
#include "mpi/common.hpp"
#include "dict.hpp"


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
    int rank = -1;
};

class LockfreeStack {
private:
    Buffer *head = nullptr;
    size_t n = 0;
public:
    void push(Buffer *bufptr)
    {
        assert(bufptr != nullptr);
        assert(bufptr->stack == nullptr);
        bufptr->stack = this;
        for (;;) {
            #pragma omp atomic read
            bufptr->stack_next = head;
            if (CAS(head, bufptr->stack_next, bufptr)) {
                #pragma omp atomic update
                n += 1;
                return;
            }
        }
    }

    Buffer * try_pop()
    {
        Buffer *bufptr;
        for (;;) {
            #pragma omp atomic read
            bufptr = head;
            if (bufptr == nullptr)
                return nullptr;
            assert(bufptr->stack == this);
            if (CAS(head, bufptr, bufptr->stack_next))
                break;
        }
        bufptr->stack = nullptr;
        #pragma omp atomic update
        n -= 1;
        return bufptr;
    }

    bool empty()
    {
        Buffer *maybe;
        #pragma omp atomic read
        maybe = head;
        return (maybe == nullptr);
    }

    size_t size()
    {
        size_t tmp;
        #pragma omp atomic read
        tmp = n;
        return tmp;
    }
};

/* Manage a collection of buffers.  The n vectors are **reserved** (but not **resized**) to capacity */
class FreeList {
private:
    LockfreeStack free;
    vector<Buffer> all;
    size_t capacity = 0;

public:
    void setup(size_t cap, size_t n)
    {
        capacity = cap;
        all.resize(n);
        for (size_t i = 0; i < n; i++) {
            all[i].reserve(capacity);
            free.push(& all[i]);
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

    void clear()
    {
        if (free.size() != all.size())
            println("freelist reset with busy buffers?");
        assert(free.size() == all.size());
    }
};



template <class Problem>
class ClawSearchProcess {
private:
    LockfreeStack outgoing;      // buffers waiting to be sent
    LockfreeStack incoming;                 // received buffers waiting to be processed
    vector<Buffer *> sendbuf, recvbuf;                // coordinator: active send/recv buffers
    vector<MPI_Request> sendreq, recvreq;           // coordinator: active send/recv requests
    bool coordinator_priority;                      // if the coordinator lacks reception buffers, the workers can't allocate

    MpiParameters &params;
    size_t capacity;                                // shorthand for params.mpi_buffer_capacity
    MPI_Comm comm;                                  // shorthand for params.comm
    u64 chunksize;                                  // shorthand for params.chunksize

    u64 N, lo, hi;
    int n_workers_started;
    int n_workers_done;
    double start;

    /*
     * coordination between hosts
     * ==========================
     * - A host is ACTIVE if it may initiate a send operation or has ongoing send operations.
     * - A host that is no longer ACTIVE is DONE.
     * - A host sends everybody else the "DONE" management message when it is DONE.
     *   (i.e., when all its send operations have completed and it won't start new ones)
     * - The whole computation is FINISHED when all hosts are DONE (all messages have been received).
     *
     * coordination between hosts
     * ==========================
     * - A worker thread is DONE when it will no longer enqueue an outgoing buffer.
     * - When all the workers in a host are DONE, and all the outgoing buffers have been completely sent, then this host is DONE.
     *
     */ 

    int n_active_hosts;
    bool signaled_termination;                      // did I tell the other hosts that I was done?
    u64 bytes_sent;

    CompactDict dict;
    
    void deactivate_host()
    {
        #pragma omp atomic update
        n_active_hosts -= 1;
    }

    void process_service_message(int source, Buffer *bufptr)
    {
        Buffer &buf = *bufptr;
        assert(buf.size() > 0);
        switch(buf[0]) {
        case DONE:
            deactivate_host();
            break;
        }
        freelist.release(bufptr);
        if (params.verbose)
            println("buffer {} released by coordinator (service)", (void *) bufptr);
    }

    void initiate_reception()
    {
        Buffer *bufptr = freelist.try_acquire();
        if (bufptr == nullptr) {
            #pragma omp atomic write
            coordinator_priority = 1;
            return;
        }
        if (params.verbose)
            println("buffer {} acquired by coordinator (for recv)", (void *) bufptr);
        #pragma omp atomic write
        coordinator_priority = 0;
        Buffer &buf = *bufptr;
        buf.resize(capacity);
        recvbuf.push_back(bufptr);
        MPI_Request &req = recvreq.emplace_back();
        MPI_Irecv(buf.data(), capacity, MPI_UINT64_T, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &req);
    }


public:
    FreeList freelist;
    const Problem &pb;
    vector<pair<u64, u64>> result;
    u64 ncoll = 0;
    int phase = 0;

    void worker_start()
    {
        #pragma omp atomic update
        n_workers_started += 1;
    }

    void worker_done()
    {
        #pragma omp atomic update
        n_workers_done += 1;
    }

    bool finished()
    {
        int tmp;
        #pragma omp atomic read
        tmp = n_active_hosts;
        return (tmp == 0) or (tmp == 1 and params.mpi_size == 1);
    }

    bool this_host_done()
    {
        int nws, nwd;
        #pragma omp atomic read
        nws = n_workers_started;
        #pragma omp atomic read
        nwd = n_workers_done;
        bool done = (nws > 0) && (nws == nwd) && outgoing.empty() && (sendbuf.empty());
        return done;
    }

    void start_phase(int i)
    {
        if (not recvbuf.empty())
            printf("warning: starting new phase with active reception buffers");
        if (not sendbuf.empty())
            printf("warning: starting new phase with active send buffers");
        freelist.clear();
        coordinator_priority = 0;

        phase = i;
        lo = params.mpi_rank * N / params.mpi_size;
        hi = (params.mpi_rank + 1) * N / params.mpi_size;
        n_active_hosts = params.mpi_size;
        n_workers_started = 0;
        n_workers_done = 0;
        signaled_termination = 0;
        start = wtime();
        bytes_sent = 0;

        std::string name[2] = {"fill", "probe"}; 
        if (params.verbose)
            println("starting phase {} ({})", i, name[i]);
    }

    ClawSearchProcess(const Problem &pb, MpiParameters &_params) : 
        params(_params), capacity(params.buffer_capacity), comm(params.comm), chunksize(params.chunksize), pb(pb)
    {
        assert(params.n_threads != 0);
        u64 n_buffers = params.send_buffers_ratio * params.mpi_size * params.n_threads + params.n_recv_buffers + params.n_slack_buffers;
        N = 1ull << pb.n;
        ncoll = 0;

        freelist.setup(capacity, n_buffers);
        dict.set_size((params.dict_capacity_ratio * N) / params.mpi_size);
        recvbuf.reserve(params.n_recv_buffers);

        if (params.verbose) {
            char hbsize[8], hdsize[8];
            human_format(n_buffers * capacity * sizeof(u64), hbsize);
            human_format(dict.nbytes(), hdsize);
            println("RAM per node == {}B buffer ({} x {}) +{}B dict ({} slots)", hbsize, n_buffers, capacity,  hdsize, dict.n_slots);
        }
    }

    void verbosity()
    {
        // worker progress ; "hash rate" ; network (MB/s, busy send, busy recv); mem (available buffers)
        u64 range_size = N / params.mpi_size;
        u64 remaining = hi - lo;
        u64 done = range_size - remaining;
        double progress = 100. * done / range_size;
        char hfrate[8], hnrate[8];
        double delta = wtime() - start;
        human_format(done / delta, hfrate);
        human_format(bytes_sent / delta, hnrate);
        println("\rDone: {:.1f}%. {} f/s. net: {}B/s ({} active send, {} ready recv, {} outgoing, {} incoming). {} free buffers", 
            progress, hfrate, hnrate, sendbuf.size(), recvbuf.size(), outgoing.size(), incoming.size(), freelist.size());
        std::fflush(stdout);
    }

    void coordinator()
    {
        // prepare the reception buffers
        for (size_t i = 0; i < recvbuf.capacity(); i++)
            initiate_reception();
            
        u64 mpi_buffer_size = params.mpi_size * (8 + MPI_BSEND_OVERHEAD);
        u8 mpi_buffer[mpi_buffer_size];
        MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);

        #pragma omp barrier            // needed otherwise the workers may steal all the buffers...

        int backoff = params.min_backoff;
        double last_verbose = start;

        while (n_active_hosts > 0) {                   // loop until the action has died down
            bool action = 0;

            // start sending pending outgoing buffers
            while (sendbuf.size() < params.max_concurrent_send) {
                Buffer *bufptr = outgoing.try_pop();
                if (bufptr == nullptr)
                    break;
                action = 1;
                sendbuf.push_back(bufptr);
                Buffer &buf = *bufptr;
                // println("coordinator, retrieving outgoing buffer of size {} for rank {}. Currently sending={}, outgoing= {}", buf.size(), rank, sendbuf.size(), outgoing.size());
                assert(not signaled_termination);
                MPI_Request &req = sendreq.emplace_back();
                MPI_Isend(buf.data(), buf.size(), MPI_UINT64_T, buf.rank, TAG_POINTS, comm, &req);
                bytes_sent += buf.size() * sizeof(u64);
            }
    
            // recycle completely sent buffers
            while (sendbuf.size() != 0) {
                int flag, index;
                MPI_Testany(sendreq.size(), sendreq.data(), &index, &flag, MPI_STATUS_IGNORE);
                if (flag == 0)
                    break;
                // spdlog::debug("coordinator, sent complete. Currently sending={}, outgoing= {}", sendbuf.size(), outgoing.size());
                action = 1;
                Buffer *bufptr = sendbuf[index];
                assert(sendreq[index] == MPI_REQUEST_NULL);
                std::swap(sendbuf[index], sendbuf.back());
                std::swap(sendreq[index], sendreq.back());
                sendbuf.pop_back();
                sendreq.pop_back();
                freelist.release(bufptr);
                if (params.verbose)
                    println("buffer {} released by coordinator (sent)", (void *) bufptr);
            }
    
            // process completely received buffers
            while (recvbuf.size() != 0) {
                int flag, index, count;
                MPI_Status status;
                MPI_Testany(recvreq.size(), recvreq.data(), &index, &flag, &status);
                if (flag == 0)
                    break;
                action = 1;
                MPI_Get_count(&status, MPI_UINT64_T, &count);
                Buffer *bufptr = recvbuf[index];
                Buffer &buf = *bufptr;
                buf.resize(count); 
                assert(recvreq[index] == MPI_REQUEST_NULL);
                std::swap(recvbuf[index], recvbuf.back());
                std::swap(recvreq[index], recvreq.back());
                recvbuf.pop_back();
                recvreq.pop_back();
                // spdlog::debug("coordinator, received a buffer.");
                if (status.MPI_TAG == TAG_MANAGEMENT) {
                    process_service_message(status.MPI_SOURCE, bufptr);
                } else {
                    incoming.push(bufptr);
                    // println("coordinator, pushed buffer to the incoming queue (size={}).", incoming.size());
                }
            }
    
            // potentially allocate new reception buffers
            for (size_t i = recvbuf.size(); i < recvbuf.capacity(); i++)
                initiate_reception();

            // detect local sender termination and notify other MPI processes
            if (this_host_done() and not signaled_termination) {
                u64 msg = DONE;
                for (int i = 0; i < params.mpi_size; i++)
                    if (i != params.mpi_rank)
                        MPI_Bsend(&msg, 1, MPI_UINT64_T, i, TAG_MANAGEMENT, params.comm);
                deactivate_host();
                signaled_termination = 1;
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

        println("rank {}, coordinator exiting the loop with {} free buffers", params.mpi_rank, freelist.size());

        // cancel the pending receives
        for (size_t i = 0; i < recvreq.size(); i++)
            MPI_Cancel(&recvreq[i]);
        MPI_Waitall(recvreq.size(), recvreq.data(), MPI_STATUSES_IGNORE);
        for (size_t i = 0; i < recvbuf.size(); i++) {
            freelist.release(recvbuf[i]);
            if (params.verbose)
                    println("buffer {} released by coordinator (after waitall)", (void *) &recvbuf[i]);
        }
        recvreq.clear();
        recvbuf.clear();

        println("rank {}, coordinator reaching the MPI barrier {} free buffers", params.mpi_rank, freelist.size());
        MPI_Barrier(params.comm);

        void *foo;
        int bar;
        MPI_Buffer_detach(&foo, &bar);
    }

    // this is specific to the direct mitm
    void process_incoming_buffer(Buffer * bufptr)
    {
        Buffer &buf = *bufptr;
        int maxkeys = 3 * pb.n;
        u64 keys[maxkeys];
        
        for (auto i = buf.begin(); i != buf.end(); i++) {
            u64 x = *i;
            i++;
            u64 z = *i;
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
                    #pragma omp atomic update
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
        freelist.release(bufptr);
        if (params.verbose)
            println("buffer {} released by process_incoming_buffer", (void *) bufptr);
    }

    bool poll_incoming()
    {
        bool action = 0;
        while (Buffer *bufptr = incoming.try_pop()) {
            process_incoming_buffer(bufptr);
            action = 1;
        }
        return action;
    }

    Buffer * get_worker_buffer()
    {
        double start = -1;
        for (;;) {
            bool action = 0;
            action = poll_incoming();
            bool must_wait;
            #pragma omp atomic read
            must_wait = coordinator_priority;
            if (not must_wait) {
                Buffer *bufptr = freelist.try_acquire();
                if (bufptr != nullptr) {
                    bufptr->resize(0);
                    // if (start >= 0)
                    //     println("worker got buffer after waiting {}s", wtime() - start);
                    if (params.verbose)
                        println("buffer {} acquired by worker", (void *) bufptr);
                    return bufptr;
                }
            }
            // backoff ?
            if (action and start >= 0) {
                // println("worker starved for {}s but got an incoming buffer", wtime() - start);
                start = -1;
            }
            if (not action and start < 0)
                start = wtime();   // start the timer
        }
    }

    pair<u64, u64> grab_chunk()
    {
        u64 local;
        #pragma omp atomic update capture
        { local = lo; lo += chunksize; }
        return pair(local, std::min(local + chunksize, hi));
    }

    void send_buffer(Buffer *bufptr, int rank)
    {
        if (rank == params.mpi_rank) {
            // fast-track: the current thread deals with it right now
            process_incoming_buffer(bufptr);
        } else {
            // slow track: the coordinator deals with it
            bufptr->rank = rank;
            outgoing.push(bufptr);
        }
    }

    // many such thread per node. This is specific to the direct mitm
    void worker()
    {
        vector<Buffer *> buffers;
        for (int i = 0; i < params.mpi_size; i++) {
            Buffer *bufptr = get_worker_buffer();
            Buffer &buf = *bufptr;
            buf.resize(0);
            buffers.push_back(bufptr);
        }

        #pragma omp barrier            // needed otherwise the workers may steal all the buffers...

        worker_start();
        for (;;) {
            poll_incoming();
            auto [chunk_lo, chunk_hi] = grab_chunk();
            if (chunk_lo >= chunk_hi)
                break;
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
            println("Running {} on MPI processes", params.mpi_size);
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
        proc.freelist.clear();

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
