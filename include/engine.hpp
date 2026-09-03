#ifndef MITM_ENGINE
#define MITM_ENGINE

#include <cassert>
#include <mpi.h>
#include <omp.h>
#include <pthread.h>
#include <err.h>
#include <memory>
#include <vector>
#include <array>

#include "counters.hpp"
#include "parameters.hpp"
#include "spsc.hpp"
#include "comm.hpp"
#include "walker.hpp"
#include "inserter.hpp"
#include "controller.hpp"

namespace mitm {

/******************************* CPU affinity ********************************/

/* pin the calling thread to `cpu`.  Returns cpu, or -1 on failure (or if cpu < 0). */
static inline int pin_to_cpu(int cpu)
{
	if (cpu < 0)
		return -1;
	cpu_set_t set;
	CPU_ZERO(&set);
	CPU_SET(cpu, &set);
	if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
		return -1;
	return cpu;
}


/******************************* the comm thread ******************************/

/*
 * One MPI rank per node, MPI_THREAD_FUNNELED.  Thread 0 is the only one that ever
 * touches MPI, and this is all of it: it drains the walker threads' queues, routes
 * points to the output buffers or straight into a local inserter's queue, services
 * the control channel, drives the end-of-round wind-down, reduces the statistics
 * and, on rank 0, runs the controller.  run() calls it at the synchronisation points
 * of a round -- begin_round(), comm_round(), end_round(), and finish() after the last
 * one -- and everything in between is here.  Its own tallies (the drops it causes)
 * go into ctx[0]->ctr like every other thread's.
 *
 * Built by thread 0 inside the OpenMP team, once pinned and once the per-thread slots
 * are published: the MPI buffers are its own (first touch, like every other thread's
 * objects), and the constructor checks the team it is about to drive.  What is shared
 * with the workers -- the parameters, the PRNG, the per-thread table and the round
 * state -- is held by reference; everything else is owned.
 */
class CommThread {
public:
	const Parameters &params;
	PRNG &prng;
	std::vector<std::unique_ptr<ThreadContext>> &ctx;     /* n_threads, indexed by tid */
	RoundState &round;

	OutBuffers outbuf;
	InBuffers inbuf;
	ControlChannel ctrl;
	Controller controller;

	bool golden_sent = false;
	bool round_over = false;           /* our TAG_END_ROUND has arrived */

	CommThread(const Parameters &params, PRNG &prng,
	           std::vector<std::unique_ptr<ThreadContext>> &ctx, RoundState &round)
		: params(params), prng(prng), ctx(ctx), round(round),
		  outbuf(params.mpi_comm, params.n_nodes, params.buffer_capacity),
		  inbuf(params.mpi_comm, params.n_in_buffers, params.buffer_capacity),
		  ctrl(params),
		  controller(params)
	{}

	/******************* routing *******************/

	/* hand a point to the inserter thread that owns its shard, on this node */
	void deliver_local(const DP &p)
	{
		int slot = (int) ((p.x / params.n_nodes) % params.inserters_per_node);
		if (not ctx[1 + slot]->q->push(p))                /* inserter `slot` is thread 1 + slot */
			ctx[0]->ctr.c[DROP_INSERTERQ] += 1;
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
					ctx[0]->ctr.c[DROP_OUT] += 1;
			}
		}
		return moved;
	}

	/*
	 * Take delivery of whatever arrived over MPI.  The Testsome/scatter/repost loop
	 * lives here rather than inside InBuffers so that the point of delivery is a
	 * direct call to deliver_local() instead of a callback.
	 */
	void poll_incoming()
	{
		int outcount = 0;
		MPI_Testsome(inbuf.req.size(), inbuf.req.data(), &outcount,
		             inbuf.done_idx.data(), inbuf.done_st.data());
		if (outcount == MPI_UNDEFINED || outcount <= 0)
			return;

		for (int t = 0; t < outcount; t++) {
			int k = inbuf.done_idx[t];
			int count = 0;
			MPI_Get_count(&inbuf.done_st[t], MPI_UINT64_T, &count);
			if (count == 0) {
				inbuf.n_sentinels += 1;       /* that node has finished the round */
			} else {
				inbuf.bytes_recv += (u64) count * sizeof(u64);
				for (int off = 0; off + DP_WORDS <= count; off += DP_WORDS) {
					DP p = {inbuf.data[k][off], inbuf.data[k][off + 1], inbuf.data[k][off + 2]};
					deliver_local(p);
				}
			}
			MPI_Irecv(inbuf.data[k].data(), inbuf.cap, MPI_UINT64_T, MPI_ANY_SOURCE,
			          TAG_POINTS, inbuf.comm, &inbuf.req[k]);
		}
	}

	/*
	 * Service the control channel: the end-of-round signal, if it came, and on rank 0
	 * the controller's inbound traffic.  This keeps running during our own drain: on
	 * rank 0 so that the reports and solutions of nodes still in their steady state are
	 * digested (round tallies, and `stop` for the next round header), and on every rank
	 * because it is the only MPI call in the last drain steps, whose progress our
	 * sentinels and rank 0's signals need to actually leave.  Nobody's liveness depends
	 * on it: every signal of the round is sent before rank 0 can enter its own drain.
	 */
	void service_control()
	{
		if (ctrl.poll()) {
			assert(not round_over);        /* exactly one signal per round */
			round_over = true;
		}
		if (params.rank == 0)
			controller.service();
	}

	/******************* one round *******************/

	/*
	 * The tallies of every thread of this rank, added up as they stand.  They are
	 * plain u64 that the workers write with no synchronisation at all: the flush
	 * makes whatever has reached memory visible, and a count that lags by a unit or
	 * two is fine for a progress report.  The exact values are read once, after the
	 * end-of-round barrier (end_round).
	 */
	void snapshot(u64 sum[N_COUNTERS])
	{
		#pragma omp flush
		for (int k = 0; k < N_COUNTERS; k++)
			sum[k] = 0;
		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N_COUNTERS; k++)
				sum[k] += ctx[t]->ctr.c[k];
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
		round.i = msg[0];
		round.root_seed = msg[1];
		round.stop = msg[2];

		if (params.rank == 0 && not round.stop)
			controller.begin_round();
	}

	/* the steady state, until our end-of-round signal arrives; then the drain */
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
			route_walker_queues();
			outbuf.poll();
			poll_incoming();

			service_control();

			/* a walker confirmed the golden collision: tell the controller now */
			if (not golden_sent && round.golden_found.load(std::memory_order_acquire)) {
				golden_sent = true;
				MPI_Bsend(round.golden, SOL_NWORDS, MPI_UINT64_T, 0, TAG_SOLUTION, params.mpi_comm);
			}

			/* periodic call home.  One-way: nothing comes back but, eventually, the
			   end-of-round signal.  Paced by the rule above. */
			if ((++poll_tick & 0xff) == 0) {
				u64 cur[N_COUNTERS];
				snapshot(cur);
				if (cur[N_DP] - prev[N_DP] < dp_report_threshold
				    && wtime() - last_ping < params.ping_delay)
					goto no_report;
				last_ping = wtime();

				/* the report is the delta of every counter since the previous one */
				u64 msg[N_COUNTERS];
				for (int k = 0; k < N_COUNTERS; k++) {
					msg[k] = cur[k] - prev[k];
					prev[k] = cur[k];
				}
				MPI_Bsend(msg, N_COUNTERS, MPI_UINT64_T, 0, TAG_REPORT, params.mpi_comm);
			}
		no_report:

			if (round_over)
				break;
		}

		drain();
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

	/*
	 * End of round.  The order matters: stop the walkers, collect everything they
	 * left behind, tell the other nodes we are done and hear the same from them, and
	 * only then let the inserters go quiet -- an inserter must never declare itself
	 * finished while points are still on their way to it.  The walkers go last,
	 * because until every inserter has stopped probing a new collision candidate can
	 * still appear.
	 */
	void drain()
	{

		/* 1. walkers stop producing points.  They keep resolving collisions. */
		set_state(WALKER, HOLD);

		/* 2. wait for them to acknowledge -- a walker only reads its state once per
		      chunk, so an empty queue before then is not an idle one -- and collect
		      whatever they left in their queues */
		for (;;) {
			bool all_held = all_in_state(WALKER, HELD);
			bool moved = route_walker_queues();
			outbuf.poll();
			poll_incoming();
			service_control();
			if (not all_held || moved)
				continue;
			bool all_empty = true;
			for (int w = 0; w < params.walkers_per_node; w++) {
				SPSCQueue &q = *ctx[1 + params.inserters_per_node + w]->q;
				if (q.head.load(std::memory_order_acquire) != q.tail.load(std::memory_order_acquire))
					all_empty = false;
			}
			if (all_empty)
				break;
		}

		/* 3. ship the partial buffers -- while still taking delivery, so the peers
		      we are waiting on can make progress too */
		while (not outbuf.flush_poll()) {
			poll_incoming();
			service_control();
		}

		/* 4. tell every node we are done sending (buffered: nothing to wait for) */
		outbuf.send_sentinels();

		/* 5. keep taking delivery until every node (ourselves included) has finished.
		      Every node was sent its end-of-round signal in one go, and rank 0 only gets
		      here after receiving its own, so nobody is waiting on us; rank 0 keeps
		      digesting reports so the round's tallies cover what the others produced
		      meanwhile. */
		while (inbuf.n_sentinels < params.n_nodes) {
			poll_incoming();
			service_control();
		}

		/* 6. nothing further will ever be delivered: the inserters can finish their
		      queues and announce themselves quiescent */
		set_state(INSERTER, DRAIN);
		while (not all_in_state(INSERTER, QUIESCENT)) {
			service_control();
			cpu_relax();
		}

		/* 7. with every inserter stopped, no new collision candidate can appear, so
		      the walkers can empty the collision queue and go quiet too.  Doing this
		      last is what stops a golden collision found on the final candidate from
		      being lost. */
		set_state(WALKER, DRAIN);
		while (not all_in_state(WALKER, QUIESCENT)) {
			service_control();
			cpu_relax();
		}
	}

	/******************* end-of-round statistics *******************/

	/*
	 * After the end-of-round barrier: every worker has returned from its round, so
	 * the tallies are exact and nobody is writing them.  Sum them over the threads
	 * of this rank, reduce over the ranks -- SUM of the counts, MAX of the
	 * HyperLogLog registers -- hand the result to the controller, and clear
	 * everything that is scoped to the round, single-threaded here by design.
	 */
	void end_round()
	{
		Counters total;
		for (int t = 0; t < params.n_threads; t++)
			total.merge(ctx[t]->ctr);

		if (params.rank == 0) {
			MPI_Reduce(MPI_IN_PLACE, total.c, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(MPI_IN_PLACE, total.hll.data(), 0x10000, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
			controller.end_round(total);
		} else {
			MPI_Reduce(total.c, NULL, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(total.hll.data(), NULL, 0x10000, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
		}

		for (int t = 0; t < params.n_threads; t++)
			ctx[t]->ctr.reset();
		inbuf.n_sentinels = 0;
	}

	/******************* termination *******************/

	/* after the last round: the answer to every rank, then close the channels */
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

		inbuf.shutdown();
		controller.shutdown();
		ctrl.shutdown();                       /* last: it detaches the Bsend buffer */
		if (not answer[0])
			return nullopt;
		return optional(tuple(answer[1], answer[2], answer[3]));
	}
};


/******************* the whole computation *******************/

/*
 * The one entry point of the engine.  Returns (i, x0, x1) -- the mixing function
 * index and the two colliding points, in wrapper coordinates -- or nothing if the
 * search gave up after opts.max_versions rounds.
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

	/* safety check: all ranks evaluate the same function */
	u64 test[3];
	u64 mask = make_mask(wrapper.m);
	test[0] = prng.rand() & mask;
	test[1] = prng.rand() & mask;
	test[2] = wrapper.mixf(test[0], test[1]);
	MPI_Bcast(test, 3, MPI_UINT64_T, 0, params.mpi_comm);
	assert(test[2] == wrapper.mixf(test[0], test[1]));

	if (params.verbose)
		Controller::banner(params, prng.seed);     /* before the shards are allocated below */

	/* What the threads share.  Each ctx slot is filled by the thread that owns it,
	   inside the team, once that thread is pinned (NUMA first touch, see there). */
	std::vector<std::unique_ptr<ThreadContext>> ctx(params.n_threads);     /* indexed by tid */
	CollisionQueue coll_q(params.coll_queue_capacity);
	RoundState round;
	optional<tuple<u64,u64,u64>> answer;
	CommThread comm(params, prng, ctx, round);

	#pragma omp parallel num_threads(params.n_threads)
	{
		int tid = omp_get_thread_num();

		/*
		 * Pin first, allocate second.  Under Linux's first-touch policy a page
		 * belongs to the NUMA node of the CPU that first writes it, so each thread
		 * builds its context -- and with it its queue and, for an inserter, its
		 * shard, whose zero-fill is the write that matters -- only once it sits on
		 * its CPU.
		 */
		int cpu = pin_to_cpu(params.thread_cpu[tid]);
		int role = WALKER;
		if (tid == 0)
			role = COMM;
		else if (tid <= params.inserters_per_node)
			role = INSERTER;
		ctx[tid] = std::make_unique<ThreadContext>(role, cpu, params);
		ThreadContext &me = *ctx[tid];

		#pragma omp barrier             // now we can inspect ctx[]

		for (;;) {
			me.state = RUNNING;
			
			#pragma omp master
			comm.begin_round(wrapper.out_mask);

			#pragma omp barrier         // makes the round header visible to all threads

			if (round.stop)
				break;

			if (me.role == COMM)
				comm.comm_round();
			else if (me.role == INSERTER)
				inserter_thread(me, params, round, coll_q);
			else
				walker_thread(me, wrapper, params, round, coll_q, tid - 1 - params.inserters_per_node);

			#pragma omp barrier

			#pragma omp master
			comm.end_round();
			
			if (me.role == INSERTER)
				me.dict->flush();

			#pragma omp barrier
		}

		#pragma omp master
		answer = comm.finish();
	}

	return answer;
}

}
#endif
