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
            if (CAS(head, bufptr, bufptr->stack_next))
                break;
        }
        assert(bufptr->stack == this);
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



template <class Problem>
class ClawSearchProcess {
private:
    FreeList intra_freelist;                       // buffers exchanged between threads
    FreeList inter_freelist;                       // buffers exchanged between nodes
    int n_active_hosts;                            // #active MPI processes
    int n_workers_started = 0;                     
    int n_workers_done = 0;
    int n_active_sends = 0;                        // #ongoing MPI send operations

    // communication between workers and coordinator
    LockfreeStack outgoing;                        // buffers waiting to be sent
    LockfreeStack incoming;                        // received buffers waiting to be processed
    
    // reception buffers
    vector<Buffer *> recvbuf;                          // active recv buffers (fixed number)
    vector<MPI_Request> recvreq;                       // corresponding requests

    MpiParameters &params;
    MPI_Comm comm;                                  // shorthand for params.comm
    u64 chunksize;                                  // shorthand for params.chunksize
    size_t bufratio;                                // shorthand for params.inter_buffer_capacity / params.inter_buffer_capacity;

    int phase = 0;
    u64 N, lo, hi;
    double start;
    u64 bytes_sent;

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
     */ 

    CompactDict dict;
    
    void process_service_message(int source, Buffer *bufptr)
    {
        Buffer &buf = *bufptr;
        assert(buf.size() > 0);
        switch(buf[0]) {
        case DONE:
            deactivate_host();
            break;
        }
    }

    void worker_start()
    {
        #pragma omp atomic update
        n_workers_started += 1;
    }

    void worker_done()
    {
        #pragma omp atomic update
        n_workers_done += 1;
        // println("mpi rank {}: {} / {} Worker done!", params.mpi_rank, n_workers_done, n_workers_started);
    }

    bool all_workers_done()
    {
        int nws, nwd;
        #pragma omp atomic read
        nws = n_workers_started;
        #pragma omp atomic read
        nwd = n_workers_done;
        return (nws > 0) && (nwd == nws);
    }

    bool coordinator_done()
    {
        int tmp;
        #pragma omp atomic read
        tmp = n_active_sends;
        return outgoing.empty() && intra_freelist.full() && (tmp == 0);
    }

    // all the send operations completed and we will no longer initiate any new send operation
    bool this_host_done()
    {
        return all_workers_done() && coordinator_done();
    }

    void deactivate_host()
    {
        #pragma omp atomic update
        n_active_hosts -= 1;
    }

    bool finished()
    {
        int tmp;
        #pragma omp atomic read
        tmp = n_active_hosts;
        return (tmp == 0) or (tmp == 1 and params.mpi_size == 1);
    }

public:
    const Problem &pb;
    vector<pair<u64, u64>> result;
    u64 ncoll = 0;


    ClawSearchProcess(const Problem &pb, MpiParameters &_params) : 
        params(_params), comm(params.comm), chunksize(params.chunksize), pb(pb)
    {
        assert(params.n_threads != 0);
        bufratio = params.inter_buffer_capacity / params.intra_buffer_capacity;
        u64 intra_inflight = params.mpi_size * params.n_threads + params.n_recv_buffers * bufratio;
        u64 n_intra_buffers = params.send_buffers_ratio * intra_inflight + params.n_slack_buffers;
        u64 n_inter_buffers = params.n_recv_buffers + params.n_slack_buffers;
        N = 1ull << pb.n;
        ncoll = 0;

        intra_freelist.setup(params.intra_buffer_capacity, n_intra_buffers);
        inter_freelist.setup(params.inter_buffer_capacity, n_inter_buffers);
        dict.set_size((params.dict_capacity_ratio * N) / params.mpi_size);
        recvbuf.reserve(params.n_recv_buffers);

        if (params.verbose) {
            char hbsize[8], hdsize[8];
            u64 intrasize = n_intra_buffers * params.intra_buffer_capacity * sizeof(u64);
            u64 intersize = n_inter_buffers * params.inter_buffer_capacity * sizeof(u64);
            human_format(intrasize + intersize, hbsize);
            human_format(dict.nbytes(), hdsize);
            println("RAM per node == {}B buffer () +{}B dict ({} slots)", hbsize, hdsize, dict.n_slots);
        }
    }


    void start_phase(int i)
    {
        intra_freelist.clear();
        inter_freelist.clear();

        phase = i;
        lo = params.mpi_rank * N / params.mpi_size;
        hi = (params.mpi_rank + 1) * N / params.mpi_size;
        n_active_hosts = params.mpi_size;
        n_workers_started = 0;
        n_workers_done = 0;
        n_active_sends = 0;
        start = wtime();
        bytes_sent = 0;

        std::string name[2] = {"fill", "probe"}; 
        if (params.verbose)
            println("starting phase {} ({})", i, name[i]);
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
        println("\rDone: {:.1f}%. {} f/s. net: {}B/s ({} active send, {} ready recv, {} incoming). {} / {} free buffers", 
            progress, hfrate, hnrate, n_active_sends, recvbuf.size(),  
            incoming.size(), intra_freelist.size(), inter_freelist.size());
        std::fflush(stdout);
    }

    void coordinator()
    {
        vector<Buffer> sendbuf;                            // send buffers (one per node)
        vector<Buffer *> recvbuf;                          // recv buffers (one per node, via the freelist)
        vector<MPI_Request> sendreq, recvreq;              // requests for the active send buffers
        vector<vector<Buffer *>> pending;                  // intra-buffers waiting for concatenation
        vector<bool> sendbusy;
        std::queue<int> empty_recv;

        bool signaled_termination = 0;                     // did I tell the other hosts that we are done?
        bool signaled_workers_done = 0;                    // did I tell the other hosts that we are done?
 
        // prepare the buffers
        sendbuf.resize(params.mpi_size);
        recvbuf.resize(params.mpi_size);
        sendreq.resize(params.mpi_size);
        recvreq.resize(params.mpi_size);
        pending.resize(params.mpi_size);
        sendbusy.resize(params.mpi_size);
        for (int i = 0; i < params.mpi_size; i++) {
            sendbuf[i].reserve(params.inter_buffer_capacity);
            recvreq[i] = MPI_REQUEST_NULL;
            sendreq[i] = MPI_REQUEST_NULL;
            sendbusy[i] = 0;
            if (i != params.mpi_rank)
                empty_recv.push(i);                // this will make the main loop allocate and initiate all the recv on the first iteration
        }

        // prepare a buffer for the MPI_Bsends;
        u64 mpi_buffer_size = params.mpi_size * (8 + MPI_BSEND_OVERHEAD);
        u8 mpi_buffer[mpi_buffer_size];
        MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);

        int backoff = params.min_backoff;
        double last_verbose = start;

        while (not finished()) {                   // loop until the action has died down
            bool action = 0;

            assert(sendreq[params.mpi_rank] == MPI_REQUEST_NULL);
            assert(recvreq[params.mpi_rank] == MPI_REQUEST_NULL);
            assert(n_active_sends >= 0);

            // recycle completely sent buffers
            for (;;) {
                int flag, i;
                MPI_Testany(sendreq.size(), sendreq.data(), &i, &flag, MPI_STATUS_IGNORE);
                if (i == MPI_UNDEFINED)
                    break;
                assert(flag != 0);
                action = 1;
                assert(not signaled_termination);
                assert(sendreq[i] == MPI_REQUEST_NULL);
                assert(sendbusy[i]);
                sendbuf[i].clear();
                sendbusy[i] = 0;
                n_active_sends -= 1;
                // println("Send to {} completed [active={}]", i, n_active_sends);
            }

            // bring intra-buffer from the workers into the pending queues
            while (Buffer *bufptr = outgoing.try_pop()) {
                assert(bufptr->rank >= 0);
                assert(bufptr->rank != params.mpi_rank);
                pending[bufptr->rank].push_back(bufptr);
                action = 1;
            }

            // process the pending queues
            bool flushing = all_workers_done();
            if (flushing and not signaled_workers_done) {
                signaled_workers_done = 1;
                prinln("MPI rank {}, all workers done\n", params.mpi_rank);
            }
            for (int i = 0; i < params.mpi_size; i++) {
                // println("to rank {}, sendbusy = {}, |pending| = {}", i, (int) sendbusy[i], pending[i].size());
                if (sendbusy[i] || pending[i].empty())
                    continue;
                if ((not flushing) and (pending[i].size() < bufratio))     // send only when large enough
                    continue;
                assert(not signaled_termination);
                action = 1;
                n_active_sends += 1;
                assert(not sendbusy[i]);
                sendbusy[i] = 1;
                for (size_t j = 0; j < bufratio; j++) {
                    if (pending[i].empty())
                        break;
                    Buffer *bufptr = pending[i].back();
                    sendbuf[i].insert(sendbuf[i].end(), bufptr->begin(), bufptr->end());
                    // for (auto x : *bufptr)
                    //     sendbuf[i].push_back(x);
                    pending[i].pop_back();
                    intra_freelist.release(bufptr);
                }
                MPI_Isend(sendbuf[i].data(), sendbuf[i].size(), MPI_UINT64_T, i, TAG_POINTS, comm, &sendreq[i]);
                // println("Send to {} started ({} items) [flushing={} / active={}]", i, sendbuf[i].size(), flushing, n_active_sends);
                bytes_sent += sendbuf[i].size() * sizeof(u64);
            }
    
            // process completely received buffers
            for (;;) {
                int flag, i, count;
                MPI_Status status;
                MPI_Testany(recvreq.size(), recvreq.data(), &i, &flag, &status);
                if (i == MPI_UNDEFINED)
                    break;
                assert(flag != 0);
                action = 1;
                MPI_Get_count(&status, MPI_UINT64_T, &count);
                // println("completed recv from rank {} ({} items)", i, count);
                Buffer *bufptr = recvbuf[i];
                assert(bufptr != nullptr);
                assert(recvreq[i] == MPI_REQUEST_NULL);
                recvbuf[i] = nullptr;
                empty_recv.push(i);
                bufptr->resize(count); 
                if (status.MPI_TAG == TAG_MANAGEMENT) {
                    process_service_message(status.MPI_SOURCE, bufptr);
                    inter_freelist.release(bufptr);
                } else {
                    incoming.push(bufptr);             // leave it to the workers
                }
            }
    
            // potentially allocate new reception buffers
            while (not empty_recv.empty()) {
                int i = empty_recv.front();
                assert(recvbuf[i] == nullptr);
                assert(recvreq[i] == MPI_REQUEST_NULL);
                Buffer *bufptr = inter_freelist.try_acquire();
                if (bufptr == nullptr)
                    break;
                action = 1;
                empty_recv.pop();
                bufptr->resize(params.inter_buffer_capacity);
                recvbuf[i] = bufptr;
                // println("starting recv from rank {}", i);
                MPI_Irecv(bufptr->data(), bufptr->size(), MPI_UINT64_T, i, MPI_ANY_TAG, comm, &recvreq[i]);
            }

            // detect local sender termination and notify other MPI processes
            if (this_host_done() and not signaled_termination) {
                u64 msg = DONE;
                // println("Signaling termination");
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

        assert(n_active_sends == 0);
        assert(outgoing.empty());
        for (int i = 0; i < params.mpi_size; i++) {
            assert(not sendbusy[i]);
            assert(pending[i].empty());
            assert(sendreq[i] == MPI_REQUEST_NULL);
        }

        // cancel the pending receives
        for (size_t i = 0; i < recvreq.size(); i++)
            if (recvreq[i] != MPI_REQUEST_NULL)
                MPI_Cancel(&recvreq[i]);
        MPI_Waitall(recvreq.size(), recvreq.data(), MPI_STATUSES_IGNORE);
        for (size_t i = 0; i < recvbuf.size(); i++)
            if (recvbuf[i] != nullptr)
                inter_freelist.release(recvbuf[i]);
        recvreq.clear();
        recvbuf.clear();
        
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
    }

    bool poll_incoming()
    {
        bool action = 0;
        while (Buffer *bufptr = incoming.try_pop()) {
            // potentially: take a bit of the work from the buffer and re-push it
            process_incoming_buffer(bufptr);
            inter_freelist.release(bufptr);
            action = 1;
        }
        return action;
    }

    Buffer * get_worker_buffer()
    {
        double start = -1;
        for (;;) {
            bool action = poll_incoming();
            Buffer *bufptr = intra_freelist.try_acquire();
            if (bufptr != nullptr) {
                bufptr->clear();
                // if (start >= 0)
                //     println("worker got buffer after waiting {}s", wtime() - start);
                // if (params.verbose)
                //     println("buffer {} acquired by worker", (void *) bufptr);
                return bufptr;
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
            intra_freelist.release(bufptr);
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
            buffers.push_back(bufptr);
        }

        worker_start();

        for (;;) {
            poll_incoming();
            auto [chunk_lo, chunk_hi] = grab_chunk();
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
