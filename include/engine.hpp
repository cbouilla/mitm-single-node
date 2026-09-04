#ifndef MITM_ENGINE
#define MITM_ENGINE

#include <cassert>
#include <mpi.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <err.h>
#include <memory>
#include <vector>
#include <array>

#include "parameters.hpp"
#include "spsc.hpp"
#include "comm.hpp"
#include "walker.hpp"
#include "inserter.hpp"
#include "controller.hpp"

namespace mitm {


/******************************* the comm thread ******************************/

/*
 * One MPI rank per node, MPI_THREAD_FUNNELED.  Thread 0 is the only one that ever
 * touches MPI, and this is all of it: it drains the walker threads' queues, routes
 * points to the output buffers or straight into a local inserter's queue, services
 * the control channel, winds the round down, reduces the statistics and, on rank 0,
 * runs the controller.  run() calls it at the synchronisation points of a round --
 * begin_round(), comm_round(), and finish() after the last one -- and everything in
 * between is here.  comm_round() is one loop, whose phase is the comm thread's own
 * ThreadContext::state (thread_state in comm.hpp), followed by end_round(), the
 * round's statistics.  Its own tallies (the drops it causes) go into ctx[0]->ctr like
 * every other thread's.
 *
 * Built by run() before the team, on the main thread, which becomes thread 0 of the
 * team: it is that thread's object, and everything it owns -- the MPI buffers, the
 * requests, the controller -- is private to it.  What it shares with the workers,
 * the parameters, the PRNG and the SharedContext, is held by reference.
 *
 * Incoming points land in a pool of n_in_buffers always-posted MPI_ANY_SOURCE
 * receives, so that buffer memory is set by that knob rather than by the number of
 * peers; poll_incoming() runs the Testsome / scatter / repost loop over the pool.
 * A zero-length TAG_POINTS message is a node's end-of-round sentinel: n_sentinels
 * counts them, and the round's incoming traffic is over when it reaches n_nodes.
 *
 * The node's end of the control channel (the three tags are listed in comm.hpp) is
 * here too, and it is small: one always-posted receive for the end-of-round signal.
 * It only ever comes from rank 0, so the receive names its source and there is no
 * wildcard at all.  The inbound side of rank 0 -- the report and solution receives --
 * belongs to the Controller, the only thing that ever reads them.  Rank 0 talks to
 * itself through MPI like any other node.  Sends on the channel are plain MPI_Bsend
 * at the call sites: the messages are short, buffered mode completes locally, so
 * there is nothing to track and the comm thread never blocks.  The buffer behind them
 * is per process and serves every Bsend of the engine -- reports, solutions, the
 * controller's end-of-round signals and the end-of-round sentinels; run() attaches it
 * before this object is built, sizes it (see there), and detaches it after finish().
 */
class CommThread {
public:
	const Parameters &params;
	PRNG &prng;
	SharedContext &shared;
	std::vector<std::unique_ptr<ThreadContext>> &ctx;     /* == shared.ctx: n_threads, indexed by tid */

	OutBuffers outbuf;

	/* the incoming bulk DP buffers: n_in_buffers receives, always posted */
	size_t in_cap;                                        /* u64 per buffer */
	std::vector<std::vector<u64>> in_data;
	std::vector<MPI_Request> in_req;
	std::vector<int> in_done_idx;                         /* MPI_Testsome output */
	std::vector<MPI_Status> in_done_st;                   /* MPI_Testsome output */
	int n_sentinels = 0;                                  /* nodes done sending this round */

	MPI_Request req_end_round = MPI_REQUEST_NULL;         /* the TAG_END_ROUND receive, always posted */
	Controller controller;

	bool golden_sent = false;
	bool round_over = false;           /* our TAG_END_ROUND has arrived */

	CommThread(const Parameters &params, PRNG &prng, SharedContext &shared)
		: params(params), prng(prng), shared(shared), ctx(shared.ctx),
		  outbuf(params.mpi_comm, params.n_nodes, params.buffer_capacity),
		  in_cap(DP_WORDS * params.buffer_capacity), in_data(params.n_in_buffers),
		  in_req(params.n_in_buffers, MPI_REQUEST_NULL),
		  in_done_idx(params.n_in_buffers), in_done_st(params.n_in_buffers),
		  controller(params)
	{
		// post receives
		for (int k = 0; k < params.n_in_buffers; k++) {
			in_data[k].resize(in_cap);
			MPI_Irecv(in_data[k].data(), in_cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS,
			          params.mpi_comm, &in_req[k]);
		}
		MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, params.mpi_comm, &req_end_round);
	}

	/******************* routing *******************/

	/* hand a point to the inserter thread that owns its shard, on this node */
	void deliver_local(const DP &p)
	{
		int slot = (int) ((p.x / params.n_nodes) % params.inserters_per_node);
		if (not ctx[1 + slot]->q->push(p))                /* inserter `slot` is thread 1 + slot */
			ctx[0]->ctr[DROP_INSERTERQ] += 1;
	}

	/* pull whatever the walker threads have produced and route it */
	bool route_walker_queues()
	{
		DP batch[64];
		bool moved = false;
		for (int s = 0; s < params.walkers_per_node; s++) {
			int t = 1 + params.inserters_per_node + s;          /* walker s is thread t */
			size_t k = ctx[t]->q->pop_bulk(batch, 64);
			if (k > 0)
				moved = true;
			for (size_t t = 0; t < k; t++) {
				int dst = (int) (batch[t].x % params.n_nodes);
				if (dst == params.rank)
					deliver_local(batch[t]);
				else if (not outbuf.push(batch[t], dst))
					ctx[0]->ctr[DROP_OUT] += 1;
			}
		}
		return moved;
	}

	/* take delivery of whatever arrived over MPI: scatter each completed buffer to the
	   local inserters, or count it as a sentinel, and repost its receive */
	void poll_incoming()
	{
		int outcount = 0;
		MPI_Testsome(in_req.size(), in_req.data(), &outcount, in_done_idx.data(), in_done_st.data());
		if (outcount == MPI_UNDEFINED || outcount <= 0)
			return;

		for (int t = 0; t < outcount; t++) {
			int k = in_done_idx[t];
			int count = 0;
			MPI_Get_count(&in_done_st[t], MPI_UINT64_T, &count);
			if (count == 0) {
				n_sentinels += 1;             /* that node has finished the round */
			} else {
				for (int off = 0; off + DP_WORDS <= count; off += DP_WORDS) {
					DP p = {in_data[k][off], in_data[k][off + 1], in_data[k][off + 2]};
					deliver_local(p);
				}
			}
			MPI_Irecv(in_data[k].data(), in_cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS,
			          params.mpi_comm, &in_req[k]);
		}
	}

	/*
	 * Service the control channel: the end-of-round signal, if it came, and on rank 0
	 * the controller's inbound traffic.  Called in every phase of the round, our own
	 * drain included: on rank 0 so that the reports and solutions of nodes still in
	 * their steady state are digested (round tallies, and `stop` for the next round
	 * header), and on every rank because its MPI_Test drives the progress our sentinels
	 * and rank 0's signals need to actually leave.  Nobody's liveness depends on it:
	 * every signal of the round is sent before rank 0 can enter its own drain.
	 */
	void service_control()
	{
		int flag = 0;
		MPI_Test(&req_end_round, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			assert(not round_over);        /* exactly one signal per round */
			round_over = true;
			MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, params.mpi_comm, &req_end_round);
		}
		if (params.rank == 0)
			controller.service();
	}

	/******************* one round *******************/

	/*
	 * The tallies of every thread of this rank, added up as they stand.  They are
	 * plain u64 that the workers write with no synchronisation at all: the flush
	 * makes whatever has reached memory visible, and a count that lags a little (one
	 * bumped once per chunk can sit in a register until the walker's next release
	 * store or mutex) is fine for a progress report.  The exact values are read once,
	 * when every worker has gone quiescent (end_round).
	 */
	void snapshot(u64 sum[N_COUNTERS])
	{
		#pragma omp flush
		for (int k = 0; k < N_COUNTERS; k++)
			sum[k] = 0;
		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N_COUNTERS; k++)
				sum[k] += ctx[t]->ctr[k];
	}

	/*
	 * The round header.  Rank 0 draws the function version and the root seed; 
	 * every rank learns them from the broadcast.
	 */
	void begin_round(u64 out_mask)
	{
		u64 msg[3];
		if (params.rank == 0) {
			msg[0] = prng.rand() & out_mask;     /* i */
			msg[1] = prng.rand();                /* root_seed */
			msg[2] = controller.stop;
		}
		MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.mpi_comm);
		shared.i = msg[0];
		shared.root_seed = msg[1];
		shared.stop = msg[2];

		if (params.rank == 0 && not shared.stop)
			controller.begin_round();
	}

	/* have all threads of `role` reached `want`? */
	bool all_in_state(int role, int want)
	{
		for (int t = 0; t < params.n_threads; t++)
			if (ctx[t]->role == role && ctx[t]->state.load(std::memory_order_acquire) != want)
				return false;
		return true;
	}

	void set_state(int role, int st)
	{
		for (int t = 0; t < params.n_threads; t++)
			if (ctx[t]->role == role)
				ctx[t]->state.store(st, std::memory_order_release);
	}

	/* is every walker queue empty?  Tested from outside, when the round winds down */
	bool walker_queues_empty()
	{
		for (int w = 0; w < params.walkers_per_node; w++) {
			SPSCQueue &q = *ctx[1 + params.inserters_per_node + w]->q;
			if (q.head.load(std::memory_order_acquire) != q.tail.load(std::memory_order_acquire))
				return false;
		}
		return true;
	}

	/*
	 * One round: the steady state, then the seven-step drain, as one loop.  Every turn
	 * does the same things whatever the phase -- route the walkers' points, take
	 * delivery over MPI, service the control channel, call home (the golden pair as
	 * soon as a walker has confirmed it, a progress report on the pacing rule) -- then
	 * tests the exit condition of the current phase, held in ctx[0]->state (see
	 * thread_state in comm.hpp), and moves on to the next one.  Nothing in the loop
	 * blocks.  Then, with every thread of the node quiescent, the round's statistics
	 * (end_round): its reductions block, and every rank reaches them.
	 *
	 * The drain phases run the full body too, and that is harmless: from FLUSHING on
	 * the walkers are held and their queues empty, so routing moves nothing; past
	 * WAITING every point message of the round has been delivered, and none of the next
	 * round can exist yet -- a node sends round r+1 points only after the round r+1
	 * broadcast returns, which needs rank 0's round r reductions to complete, i.e.
	 * every node past its own drain.
	 */
	void comm_round()
	{
		double last_ping = wtime();
		round_over = false;

		u64 prev[N_COUNTERS] = {0};        /* the tallies at the previous report */

		/*
		 * Reporting on a timer alone lets a short round run far past beta*w: the
		 * controller cannot end it before the first reports land, and by then the
		 * walkers have produced a whole ping_delay of points.  That overshoot also
		 * runs the chain counter j past the width the dictionary encodes it in.  So
		 * also report once a node has produced its share of the round.
		 */
		u64 dp_report_threshold = params.points_per_version
		                        / ((u64) params.reports_per_round * params.n_nodes);
		if (dp_report_threshold < 1)
			dp_report_threshold = 1;
		u64 poll_tick = 0;

		for (;;) {
			bool moved = route_walker_queues();
			poll_incoming();
			service_control();

			/* a walker confirmed the golden collision: tell the controller now */
			if (not golden_sent && shared.golden_found.load(std::memory_order_acquire)) {
				MPI_Bsend(shared.golden, SOL_NWORDS, MPI_UINT64_T, 0, TAG_SOLUTION, params.mpi_comm);
				golden_sent = true;
			}

			/* periodic call home */
			if ((++poll_tick & 0xff) == 0) {
				u64 cur[N_COUNTERS];
				snapshot(cur);
				bool due = cur[N_DP] - prev[N_DP] >= dp_report_threshold || wtime() - last_ping >= params.ping_delay;
				if (due) {
					last_ping = wtime();
					/* the report is the delta of every counter since the previous one */
					u64 msg[N_COUNTERS];
					for (int k = 0; k < N_COUNTERS; k++) {
						msg[k] = cur[k] - prev[k];
						prev[k] = cur[k];
					}
					MPI_Bsend(msg, N_COUNTERS, MPI_UINT64_T, 0, TAG_REPORT, params.mpi_comm);
				}
			}

			/* the phase: its exit test, and the step it triggers.  Ours alone, hence relaxed */
			int st = ctx[0]->state.load(std::memory_order_relaxed);
			assert(st == RUNNING || st == COLLECTING || not moved);      /* held walkers produce nothing */
			
			switch (st) {
			case RUNNING:
				if (round_over) {
				/* 1. our end-of-round signal: the walkers stop producing points.  
				      They keep resolving collisions. */
					set_state(WALKER, HOLD);
					st = COLLECTING;
				}
				break;

			case COLLECTING:
				/* 2. every walker has acknowledged -- one only reads its state once per
				      chunk, so an empty queue before then is not an idle one -- and their
				      queues are empty: every point of ours is routed.  A walker's last
				      push happens-before its HELD store, so the queue test that follows
				      the acquire load of HELD sees it. */
				if (all_in_state(WALKER, HELD) && not moved && walker_queues_empty())
					st = FLUSHING;
				break;

			case FLUSHING:
				/* 3. the partial buffers are out and every send has completed -- while
				      still taking delivery, so the peers we wait on make progress too.
				   4. tell every node we are done sending (buffered: nothing to wait for) */
				if (outbuf.flush_poll()) {
					outbuf.send_sentinels();
					st = WAITING;
				}
				break;

			case WAITING:
				/* 5. every node (ourselves included) has finished sending to us.  Every
				      node was sent its end-of-round signal in one go, and rank 0 only
				      gets here after receiving its own, so nobody is waiting on us; rank 0
				      keeps digesting reports so the round's tallies cover what the others
				      produced meanwhile.
				   6. nothing further will ever be delivered: the inserters can finish
				      their queues and announce themselves quiescent */
				if (n_sentinels == params.n_nodes) {
					set_state(INSERTER, DRAIN);
					st = DRAINING_INSERTERS;
				}
				break;

			case DRAINING_INSERTERS:
				/* 6. every inserter has: every point of the round has been probed.
				   7. with every inserter stopped, no new collision candidate can appear,
				      so the walkers can empty the collision queue and go quiet too.
				      Doing this last is what stops a golden collision found on the final
				      candidate from being lost. */
				if (all_in_state(INSERTER, QUIESCENT)) {
					set_state(WALKER, DRAIN);
					st = DRAINING_WALKERS;
				}
				break;

			case DRAINING_WALKERS:
				/* 7. every walker has: every candidate of the round has been resolved,
				      and the round is over for this node */
				if (all_in_state(WALKER, QUIESCENT))
					st = QUIESCENT;
				break;
			}
			ctx[0]->state.store(st, std::memory_order_relaxed);
			if (st == QUIESCENT)
				break;
		}
		end_round();
	}

	/******************* end-of-round statistics *******************/

	/*
	 * The round's statistics, once every worker of the node is QUIESCENT: their
	 * tallies are then exact and nobody writes them, and no walker touches the
	 * HyperLogLog any more -- a worker's last increment happens-before its QUIESCENT
	 * store, which comm_round() loaded with acquire.  Sum the tallies over the threads
	 * of this rank, copy the node's registers out, reduce both over the ranks -- SUM of
	 * the counts, MAX of the registers -- hand the result to the controller, and clear
	 * everything that is scoped to the round; nobody touches any of it again before
	 * the next round's barrier.  The reductions block, and they are the first thing a
	 * rank does once its drain is over, so every rank reaches them.
	 */
	void end_round()
	{
		u64 sum[N_COUNTERS];
		snapshot(sum);
		u8 hll[HLL_REGISTERS];
		for (int k = 0; k < HLL_REGISTERS; k++)
			hll[k] = shared.hll[k].load(std::memory_order_relaxed);

		if (params.rank == 0) {
			MPI_Reduce(MPI_IN_PLACE, sum, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(MPI_IN_PLACE, hll, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
			controller.end_round(sum, hll);
		} else {
			MPI_Reduce(sum, NULL, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(hll, NULL, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
		}

		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N_COUNTERS; k++)
				ctx[t]->ctr[k] = 0;
		for (int k = 0; k < HLL_REGISTERS; k++)
			shared.hll[k].store(0, std::memory_order_relaxed);
		n_sentinels = 0;
	}

	/******************* termination *******************/

	/* after the last round: the answer to every rank, then cancel every posted receive */
	optional<tuple<u64,u64,u64>> finish()
	{
		/* answer[0] says whether the other three mean anything */
		u64 answer[4] = {0, 0, 0, 0};
		if (params.rank == 0) {
			if (controller.solution) {
				auto [i, x0, x1] = *controller.solution;
				answer[0] = 1; answer[1] = i; answer[2] = x0; answer[3] = x1;
			}
			controller.done();
		}
		MPI_Bcast(answer, 4, MPI_UINT64_T, 0, params.mpi_comm);

		for (size_t k = 0; k < in_req.size(); k++)
			if (in_req[k] != MPI_REQUEST_NULL) {
				MPI_Cancel(&in_req[k]);
				MPI_Wait(&in_req[k], MPI_STATUS_IGNORE);
			}
		controller.shutdown();
		if (req_end_round != MPI_REQUEST_NULL) {
			MPI_Cancel(&req_end_round);
			MPI_Wait(&req_end_round, MPI_STATUS_IGNORE);
		}
		if (not answer[0])
			return nullopt;
		return optional(tuple(answer[1], answer[2], answer[3]));
	}
};


/******************* the whole computation *******************/

/*
 * The one entry point of the engine.  Returns (i, x0, x1) --- the mixing
 * function index and the two colliding points --- or nothing if the search
 * gave up after opts.max_versions rounds.
 */
template <class ProblemWrapper>
optional<tuple<u64,u64,u64>> run(const ProblemWrapper &wrapper, u64 nbytes_memory,
                                 const Options &opts, PRNG &prng)
{
	/* MPI must be initialised with at least FUNNELED support */
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI: MPI_THREAD_FUNNELED is required (got %d).  Use MPI_Init_thread.", provided);

	const Parameters params(opts, nbytes_memory, wrapper.n, wrapper.m);

	/* The MPI_Bsend buffer of the engine (preserving the previous one) */
	static_assert((int) N_COUNTERS >= (int) SOL_NWORDS, "the Bsend slots are sized for a report");
	void *prev_bsend_buf = NULL;
	int prev_bsend_size = 0;
	MPI_Buffer_detach(&prev_bsend_buf, &prev_bsend_size);
	size_t slots = 2 * (size_t) params.n_nodes + params.bsend_slack;
	size_t msg = N_COUNTERS * sizeof(u64) + MPI_BSEND_OVERHEAD;
	std::vector<char> bsend_buf(slots * msg);
	MPI_Buffer_attach(bsend_buf.data(), (int) bsend_buf.size());

	/* safety check: all ranks evaluate the same function (no uninitialized state) */
	u64 test[3];
	u64 mask = make_mask(wrapper.m);
	test[0] = prng.rand() & mask;
	test[1] = prng.rand() & mask;
	test[2] = wrapper.mixf(test[0], test[1]);
	MPI_Bcast(test, 3, MPI_UINT64_T, 0, params.mpi_comm);
	assert(test[2] == wrapper.mixf(test[0], test[1]));

	if (params.verbose)
		Controller::banner(params, prng.seed);

	SharedContext shared(params);
	CommThread comm(params, prng, shared);
	int n_unpinned = 0;                 /* threads not where Parameters put them: fatal */

	#pragma omp parallel num_threads(params.n_threads)
	{
		int tid = omp_get_thread_num();

		/*
		 * Pin first, allocate second.  Under Linux's first-touch policy a page
		 * belongs to the NUMA node of the CPU that first writes it.  So ask the
		 * kernel where we are before touching anything: a thread that is not on
		 * its CPU would void the placement, and the whole rank stops instead --
		 * through thread 0, the only one allowed to call MPI (FUNNELED), and
		 * before anyone zero-fills a shard.
		 */
		int want = params.thread_cpu[tid];
		if (want >= 0 && pin_to_cpu(want) < 0) {
			warn("MPI: rank %d: cannot pin thread %d to CPU %d", params.rank, tid, want);
			#pragma omp atomic
			n_unpinned += 1;
		}
		unsigned cpu = 0, numa_node = 0;
		if (syscall(SYS_getcpu, &cpu, &numa_node, NULL) != 0)
			err(1, "getcpu");
		if (want >= 0 && (int) cpu != want) {
			warnx("MPI: rank %d: thread %d asked for CPU %d, runs on CPU %u", params.rank, tid, want, cpu);
			#pragma omp atomic
			n_unpinned += 1;
		}

		#pragma omp barrier             /* everybody is pinned, or nobody allocates */

		if (tid == 0 && n_unpinned > 0)
			MPI_Abort(params.mpi_comm, 1);
		int role = WALKER;
		if (tid == 0)
			role = COMM;
		else if (tid <= params.inserters_per_node)
			role = INSERTER;
		shared.ctx[tid] = std::make_unique<ThreadContext>(role, cpu, numa_node, params);
		if (role == INSERTER)
			shared.shards[tid - 1] = std::make_unique<PcsDict>(params.jbits, params.w_shard);
		ThreadContext &me = *shared.ctx[tid];

		#pragma omp barrier             /* now we can inspect the shared context */

		if (tid == 0 && params.verbose)
			Controller::placement(params, shared);

		for (;;) {
			me.state = RUNNING;
			
			#pragma omp master
			comm.begin_round(wrapper.out_mask);

			#pragma omp barrier         /* makes the round header visible to all threads */

			if (shared.stop)
				break;

			if (me.role == COMM)
				comm.comm_round();          /* ends with end_round(), the statistics */
			else if (me.role == INSERTER)
				inserter_thread(me, params, shared, tid - 1);
			else
				walker_thread(me, wrapper, params, shared, tid - 1 - params.inserters_per_node);

			#pragma omp barrier         /* every thread is back.  Nothing depends on it */

			if (me.role == INSERTER)        // later: collective flush
				shared.shards[tid - 1]->flush();
		}
	}

	optional<tuple<u64,u64,u64>> answer = comm.finish();

	/* detach waits for our last buffered sends to be out; then give MPI back the
	   caller's buffer, if there was one */
	void *ptr = NULL;
	int sz = 0;
	MPI_Buffer_detach(&ptr, &sz);
	assert(ptr == bsend_buf.data());
	if (prev_bsend_size > 0)
		MPI_Buffer_attach(prev_bsend_buf, prev_bsend_size);
	return answer;
}

}
#endif
