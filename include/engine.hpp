#ifndef MITM_ENGINE
#define MITM_ENGINE

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


/********************************** the node *********************************/

/*
 * One MPI rank per node, MPI_THREAD_FUNNELED.  Thread 0 is the only one that ever
 * touches MPI: it drains the walker threads' queues, routes points to output buffers
 * or straight into a local inserter's queue, and services the control channel.  On
 * rank 0 it additionally runs the controller.
 */
template <class ProblemWrapper>
class PcsNode {
public:
	const Parameters &params;
	PRNG &prng;
	const ProblemWrapper &wrapper;                         /* stateless: every walker shares it */

	std::vector<std::unique_ptr<ThreadContext>> ctx;       /* n_threads, indexed by tid */
	std::vector<std::unique_ptr<SPSCQueue>> walker_q;    /* one per walker thread */
	std::vector<std::unique_ptr<SPSCQueue>> inserter_q;    /* one per inserter thread */
	std::vector<std::unique_ptr<PcsDict>> dict;            /* one shard per inserter thread */

	CollisionQueue coll_q;
	RoundState round;
	OutBuffers outbuf;
	InBuffers inbuf;
	ControlChannel ctrl;
	Controller controller;

	u64 n_drop_inserterq = 0;          /* comm thread could not hand a DP to a inserter */
	bool golden_sent = false;

	/* control-loop state, shared between comm_round() and drain() */
	bool awaiting_assignment = false;
	bool got_new_version = false;

	PcsNode(const ProblemWrapper &wrapper, const Parameters &params, PRNG &prng)
		: params(params), prng(prng), wrapper(wrapper),
		  coll_q(params.coll_queue_capacity),
		  outbuf(params.mpi_comm, params.n_nodes, params.buffer_capacity),
		  inbuf(params.mpi_comm, params.n_in_buffers, params.buffer_capacity),
		  ctrl(params),
		  controller(params)
	{
		if (params.verbose)
			controller.banner(prng.seed);      /* before anything big is allocated */

		for (int t = 0; t < params.n_threads; t++)
			ctx.push_back(std::make_unique<ThreadContext>());

		for (int r = 0; r < params.inserters_per_node; r++) {
			dict.push_back(std::make_unique<PcsDict>(params.jbits, params.w_shard));
			inserter_q.push_back(std::make_unique<SPSCQueue>(params.inserter_queue_capacity));
			ctx[1 + r]->role = INSERTER;
		}
		for (int w = 0; w < params.walkers_per_node; w++) {
			walker_q.push_back(std::make_unique<SPSCQueue>(params.walker_queue_capacity));
			ctx[1 + params.inserters_per_node + w]->role = WALKER;
		}
		ctx[0]->role = COMM;
	}


	/* hand a point to the inserter thread that owns its shard, on this node */
	void deliver_local(const DP &p)
	{
		int slot = (int) ((p.x / params.n_nodes) % params.inserters_per_node);
		if (not inserter_q[slot]->push(p))
			n_drop_inserterq += 1;
	}

	/* pull whatever the walker threads have produced and route it */
	bool route_walker_queues()
	{
		DP batch[64];
		bool moved = false;
		for (int s = 0; s < params.walkers_per_node; s++) {
			size_t k = walker_q[s]->pop_bulk(batch, 64);
			if (k > 0)
				moved = true;
			for (size_t t = 0; t < k; t++) {
				int dst = (int) (batch[t].x % params.n_nodes);
				if (dst == params.rank)
					deliver_local(batch[t]);
				else
					outbuf.push(batch[t], dst);       /* counts its own drops */
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
	 * Service the control channel.  This has to keep running during our own drain:
	 * rank 0 hands itself NEW_VERSION like any other node and would otherwise stop
	 * answering everyone else, who would then never stop producing.
	 */
	void service_control()
	{
		u64 msg[REP_NWORDS];
		int src = -1;
		int count = ctrl.poll(msg, &src);

		if (count == 0)
			return;

		if (count == 1) {                       /* an assignment, for us */
			awaiting_assignment = false;
			if (msg[0] == NEW_VERSION)
				got_new_version = true;
			return;
		}

		/* a node report; only the controller is ever sent one */
		int a = controller.handle_report(msg);
		if (a >= 0) {
			u64 reply = (u64) a;
			MPI_Bsend(&reply, 1, MPI_UINT64_T, src, TAG_CONTROL, params.mpi_comm);
		}
	}

	/******************* the comm thread, for one round *******************/

	void comm_round()
	{

		double last_ping = wtime();
		awaiting_assignment = false;
		got_new_version = false;

		u64 prev_dp = 0, prev_ds = 0, prev_dc = 0, prev_do = 0, prev_dr = 0, prev_np = 0;

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
				u64 msg[REP_NWORDS] = {0};
				msg[REP_GOLDEN] = 1;
				msg[REP_I] = round.golden[0];
				msg[REP_X0] = round.golden[1];
				msg[REP_X1] = round.golden[2];
				MPI_Bsend(msg, REP_NWORDS, MPI_UINT64_T, 0, TAG_CONTROL, params.mpi_comm);
			}

			/* periodic call home.  Only ever one outstanding, so the controller can
			   count each node exactly once when it switches to a new version. */
			if (not awaiting_assignment && ((++poll_tick & 0xff) == 0)) {
				u64 cur_dp = 0, cur_ds = 0, cur_dc = 0, cur_np = 0;
				for (int s = 0; s < params.walkers_per_node; s++) {
					int t = 1 + params.inserters_per_node + s;
					cur_dp += ctx[t]->n_dp.load(std::memory_order_relaxed);
					cur_ds += ctx[t]->n_drop_walkerq.load(std::memory_order_relaxed);
				}
				for (int r = 0; r < params.inserters_per_node; r++) {
					cur_dc += ctx[1 + r]->n_drop_coll.load(std::memory_order_relaxed);
					cur_np += ctx[1 + r]->n_probe.load(std::memory_order_relaxed);
				}
				u64 cur_do = outbuf.n_dropped;
				u64 cur_dr = n_drop_inserterq;

				if (cur_dp - prev_dp < dp_report_threshold
				    && wtime() - last_ping < params.ping_delay)
					goto no_report;
				last_ping = wtime();

				u64 msg[REP_NWORDS] = {0};
				msg[REP_NDP] = cur_dp - prev_dp;
				msg[REP_DROP_WALKERQ] = cur_ds - prev_ds;
				msg[REP_DROP_OUT] = cur_do - prev_do;
				msg[REP_DROP_INSERTERQ] = cur_dr - prev_dr;
				msg[REP_DROP_COLL] = cur_dc - prev_dc;
				msg[REP_NPROBE] = cur_np - prev_np;
				prev_dp = cur_dp; prev_ds = cur_ds; prev_dc = cur_dc;
				prev_do = cur_do; prev_dr = cur_dr; prev_np = cur_np;

				MPI_Bsend(msg, REP_NWORDS, MPI_UINT64_T, 0, TAG_CONTROL, params.mpi_comm);
				awaiting_assignment = true;
			}
		no_report:

			if (got_new_version)
				break;
		}

		drain();
	}

	/*
	 * End of round.  The order here matters: stop the walkers, empty everything they
	 * left behind, tell the other nodes we are done, keep receiving until they all say
	 * the same, then let the inserters and the collision queue run dry.
	 */
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
			for (int w = 0; w < params.walkers_per_node; w++)
				if (walker_q[w]->head.load(std::memory_order_acquire)
				    != walker_q[w]->tail.load(std::memory_order_acquire))
					all_empty = false;
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
		      The other nodes only get here once the controller has told them to, so
		      rank 0 must keep answering reports throughout. */
		while (inbuf.n_sentinels < params.n_nodes) {
			poll_incoming();
			service_control();
		}

		/* 5b. rank 0 must not leave the round until every node has been handed its
		       NEW_VERSION -- a node still waiting for a reply would never reach the
		       collective epilogue, and the whole job would deadlock there. */
		if (params.rank == 0)
			while (controller.n_active_nodes > 0) {
				service_control();
				cpu_relax();
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

	void epilogue()
	{
		Counters total;
		for (int t = 0; t < params.n_threads; t++)
			total.merge(ctx[t]->ctr);

		/* everything in here is scoped to the round: reset() clears it below */
		u64 st[ST_NWORDS];
		st[ST_NEVAL] = total.n_eval;
		st[ST_NPOINTS_TRAILS] = total.n_points_trails;
		st[ST_NCOLL] = total.n_collisions;
		st[ST_LEN_MIN] = total.colliding_len_min;
		st[ST_LEN_MAX] = total.colliding_len_max;
		st[ST_BAD_PROBE] = total.bad_probe;
		st[ST_BAD_ROBINHOOD] = total.bad_walk_robinhood;
		st[ST_BAD_NONCOLLIDING] = total.bad_walk_noncolliding;
		st[ST_BAD_COLLISION] = total.bad_collision;
		st[ST_BAD_DP] = total.bad_dp;

		if (params.rank == 0)
			MPI_Reduce(MPI_IN_PLACE, st, ST_NWORDS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
		else
			MPI_Reduce(st, NULL, ST_NWORDS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);

		vector<u8> hll = total.hll;
		if (params.rank == 0)
			MPI_Reduce(MPI_IN_PLACE, hll.data(), 0x10000, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
		else
			MPI_Reduce(hll.data(), NULL, 0x10000, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);

		if (params.rank == 0)
			controller.end_round(st, hll);

		/* reset every thread's per-round state (single-threaded here, by design) */
		for (int t = 0; t < params.n_threads; t++) {
			ctx[t]->ctr.reset();
			ctx[t]->n_dp.store(0, std::memory_order_relaxed);
			ctx[t]->n_probe.store(0, std::memory_order_relaxed);
			ctx[t]->n_drop_walkerq.store(0, std::memory_order_relaxed);
			ctx[t]->n_drop_coll.store(0, std::memory_order_relaxed);
		}
		outbuf.n_dropped = 0;
		n_drop_inserterq = 0;
		inbuf.n_sentinels = 0;
	}

	/******************* the whole computation *******************/

	optional<tuple<u64,u64,u64>> run()
	{
		#pragma omp parallel num_threads(params.n_threads)
		{
			int tid = omp_get_thread_num();
			ThreadContext &me = *ctx[tid];
			me.cpu = pin_to_cpu(params.thread_cpu[tid]);

			#pragma omp barrier
			#pragma omp master
			if (params.verbose && params.bind_threads) {
				printf("MPI: rank 0 thread->cpu map:");
				for (int t = 0; t < params.n_threads; t++)
					printf(" %d:%s%d", t,
						(t == 0) ? "c" : (t <= params.inserters_per_node ? "i" : "w"),
						ctx[t]->cpu);
				printf("\n");
				fflush(stdout);
			}

			for (;;) {
				me.state = RUNNING;

				#pragma omp master
				{
					// acquire the new function version from the controller
					u64 msg[3];
					if (params.rank == 0) {
						msg[0] = prng.rand() & wrapper.out_mask;     // i
						msg[1] = prng.rand();                // root_seed
						msg[2] = controller.stop;
					}
					MPI_Bcast(msg, 3, MPI_UINT64_T, 0, params.mpi_comm);
					round.i = msg[0];
					round.root_seed = msg[1];
					round.stop = msg[2];						
					/* no point arming a round that nobody will run: everyone breaks
					   out just past the barrier below */
					if (params.rank == 0 && not round.stop)
						controller.begin_round();
				}

				#pragma omp barrier         // makes the result of the Bcast visible to all threads

				if (round.stop)
					break;

				if (tid == 0)
					comm_round();
				else if (me.role == INSERTER)
					inserter_thread(me, params, round, *dict[tid - 1],
					                *inserter_q[tid - 1], coll_q);
				else
					walker_thread(me, wrapper, params, round,
					              *walker_q[tid - 1 - params.inserters_per_node], coll_q,
					              tid - 1 - params.inserters_per_node);

				#pragma omp barrier

				if (tid == 0)
					epilogue();
				else if (me.role == INSERTER)
					dict[tid - 1]->flush();

				#pragma omp barrier
			}
		}

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
		ctrl.shutdown();
		if (not answer[0])
			return nullopt;
		return optional(tuple(answer[1], answer[2], answer[3]));
	}
};


/*
 * The one entry point of the engine.  Returns (i, x0, x1) -- the mixing function
 * index and the two colliding points, in wrapper coordinates -- or nothing if the
 * search gave up after opts.max_versions rounds.
 *
 * This is where the Options, the RAM budget and the wrapped problem's dimensions
 * meet, so this is where the Parameters are derived (once, and for everything below).
 */
template<class ProblemWrapper>
optional<tuple<u64,u64,u64>> run_engine(const ProblemWrapper &wrapper, u64 nbytes_memory,
                                        const Options &opts, PRNG &prng)
{
	/* MPI must be initialised with at least FUNNELED support */
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI: MPI_THREAD_FUNNELED is required (got %d).  Use MPI_Init_thread.", provided);

	Parameters params(opts, nbytes_memory, wrapper.n, wrapper.m);

	/* safety check: all ranks evaluate the same function */
	u64 test[3];
	u64 mask = make_mask(wrapper.m);
	test[0] = prng.rand() & mask;
	test[1] = prng.rand() & mask;
	test[2] = wrapper.mixf(test[0], test[1]);
	MPI_Bcast(test, 3, MPI_UINT64_T, 0, params.mpi_comm);
	assert(test[2] == wrapper.mixf(test[0], test[1]));

	PcsNode<ProblemWrapper> node(wrapper, params, prng);
	return node.run();
}

}
#endif
