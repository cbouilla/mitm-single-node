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
 * Thread 0's object: all of the rank's MPI.  Routing, the receive pool, the control channel, the
 * round's state machine (comm_round), the end-of-round reductions and, on rank 0, the controller.
 * Built by run() on the main thread, before the team.  PROTOCOL.md §2 and §4.
 */
class CommThread {
public:
	const Parameters &params;
	PRNG &prng;
	SharedContext &shared;
	std::vector<std::unique_ptr<ThreadContext>> &ctx;     /* == shared.ctx: n_threads, indexed by tid */

	OutBuffers outbuf;

	/* the receive pool: n_in_buffers always-posted ANY_SOURCE receives, sized by that knob and not by the peers */
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
		bucket.resize((size_t) params.dicts_per_node * BUCKET_CAP);
		bucket_n.assign(params.dicts_per_node, 0);
		for (int k = 0; k < params.n_in_buffers; k++) {
			in_data[k].resize(in_cap);
			MPI_Irecv(in_data[k].data(), in_cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS,
			          params.mpi_comm, &in_req[k]);
		}
		MPI_Irecv(NULL, 0, MPI_UINT64_T, 0, TAG_END_ROUND, params.mpi_comm, &req_end_round);
	}

	/******************* routing *******************/

	/*
	 * A point is staged in its destination shard's bucket and handed to the dict thread in runs, never
	 * one at a time: a push per point is one release store on a line the consumer is spinning on, and
	 * one sweep of the producer queues scatters over every shard of the node, so each ring's lines
	 * ping-pong between the two cores once per point (PROBLEM.md §3, §5.C).  A bucket never outlives
	 * the call that filled it -- route_producer_queues() and poll_incoming() both flush before
	 * returning -- so the drain of PROTOCOL.md §4.5 sees exactly what it did before.
	 */
	static constexpr size_t BUCKET_CAP = 128;             /* points staged per shard */
	std::vector<DP> bucket;                               /* dicts_per_node runs of BUCKET_CAP */
	std::vector<size_t> bucket_n;                         /* points staged for each shard */

	/* hand shard `slot`'s staged run to its dict thread, dropping whatever did not fit */
	void flush_bucket(int slot)
	{
		size_t n = bucket_n[slot];
		if (n == 0)
			return;
		size_t k = ctx[1 + slot]->q->push_bulk(&bucket[slot * BUCKET_CAP], n);   /* dict thread `slot` is thread 1 + slot */
		ctx[0]->ctr[DROP_DICTQ] += n - k;
		bucket_n[slot] = 0;
	}

	void flush_local()
	{
		for (int slot = 0; slot < params.dicts_per_node; slot++)
			flush_bucket(slot);
	}

	/* stage a point for the dict thread that owns its shard, on this node */
	void deliver_local(const DP &p)
	{
		int slot = (int) ((p.x / params.n_nodes) % params.dicts_per_node);
		if (bucket_n[slot] == BUCKET_CAP)
			flush_bucket(slot);
		bucket[slot * BUCKET_CAP + bucket_n[slot]] = p;
		bucket_n[slot] += 1;
	}

	/* pull whatever the producer threads have produced and route it */
	bool route_producer_queues()
	{
		DP batch[64];
		bool moved = false;
		for (int s = 0; s < params.producers_per_node; s++) {
			int t = 1 + params.dicts_per_node + s;              /* producer s is thread t */
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
		flush_local();
		return moved;
	}

	/* take delivery of whatever arrived over MPI, and repost each completed receive */
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
					DP p = {in_data[k][off], in_data[k][off + 1]};
					deliver_local(p);
				}
			}
			MPI_Irecv(in_data[k].data(), in_cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS,
			          params.mpi_comm, &in_req[k]);
		}
		flush_local();
	}

	/*
	 * The control channel: our end-of-round signal and, on rank 0, the controller's receives.  Run in
	 * every phase of the round, our own drain included: PROTOCOL.md §4.5, last paragraph.
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

	/* the node's tallies as they stand: approximate mid-round, exact once every worker is QUIESCENT (§3.4) */
	void snapshot(u64 sum[N_COUNTERS])
	{
		#pragma omp flush
		for (int k = 0; k < N_COUNTERS; k++)
			sum[k] = 0;
		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N_COUNTERS; k++)
				sum[k] += ctx[t]->ctr[k];
	}

	/* the round header: drawn by rank 0, broadcast, stored in the SharedContext (PROTOCOL.md §4.2) */
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

	/* ask every thread of `role` to move to `st` */
	void set_state(int role, int st)
	{
		for (int t = 0; t < params.n_threads; t++)
			if (ctx[t]->role == role)
				ctx[t]->state.store(st, std::memory_order_release);
	}

	/* is every producer queue empty?  Tested from outside, when the round winds down */
	bool producer_queues_empty()
	{
		for (int w = 0; w < params.producers_per_node; w++) {
			SPSCQueue &q = *ctx[1 + params.dicts_per_node + w]->q;
			if (q.head.load(std::memory_order_acquire) != q.tail.load(std::memory_order_acquire))
				return false;
		}
		return true;
	}

	/*
	 * One round of the comm thread: steady state and drain as one non-blocking loop, whose phase is
	 * ctx[0]->state, then end_round().  The turn's body and each phase's exit test: PROTOCOL.md §4.3, §4.5.
	 */
	void comm_round()
	{
		double last_ping = wtime();
		round_over = false;

		u64 prev[N_COUNTERS] = {0};        /* the tallies at the previous report */

		/* report by volume too: on a timer alone a short round overshoots beta*w and runs j past jbits */
		u64 dp_report_threshold = params.points_per_version
		                        / ((u64) params.reports_per_round * params.n_nodes);
		if (dp_report_threshold < 1)
			dp_report_threshold = 1;
		u64 poll_tick = 0;

		for (;;) {
			bool moved = route_producer_queues();
			poll_incoming();
			service_control();

			if (not golden_sent && shared.golden_found.load(std::memory_order_acquire)) {
				MPI_Bsend(shared.golden, SOL_NWORDS, MPI_UINT64_T, 0, TAG_SOLUTION, params.mpi_comm);
				golden_sent = true;
			}

			if ((++poll_tick & 0xff) == 0) {
				u64 cur[N_COUNTERS];
				snapshot(cur);
				bool due = cur[N_DP] - prev[N_DP] >= dp_report_threshold || wtime() - last_ping >= params.ping_delay;
				if (due) {
					last_ping = wtime();
					u64 msg[N_COUNTERS];
					for (int k = 0; k < N_COUNTERS; k++) {
						msg[k] = cur[k] - prev[k];
						prev[k] = cur[k];
					}
					MPI_Bsend(msg, N_COUNTERS, MPI_UINT64_T, 0, TAG_REPORT, params.mpi_comm);
				}
			}

			/* our phase: nobody else reads it, hence relaxed */
			int st = ctx[0]->state.load(std::memory_order_relaxed);
			assert(st == RUNNING || st == COLLECTING || not moved);      /* held producers produce nothing */

			switch (st) {
			case RUNNING:
				if (round_over) {                                /* step 1 */
					set_state(PRODUCER, HOLD);
					st = COLLECTING;
				}
				break;

			case COLLECTING:                                     /* step 2 */
				if (all_in_state(PRODUCER, HELD) && not moved && producer_queues_empty())
					st = FLUSHING;
				break;

			case FLUSHING:                                       /* steps 3, 4 */
				if (outbuf.flush_poll()) {
					outbuf.send_sentinels();
					st = WAITING;
				}
				break;

			case WAITING:                                        /* steps 5, 6 */
				if (n_sentinels == params.n_nodes) {
					set_state(DICT, DRAIN);
					st = DRAINING_DICTS;
				}
				break;

			case DRAINING_DICTS:                                 /* steps 6, 7 */
				if (all_in_state(DICT, QUIESCENT)) {
					set_state(PRODUCER, DRAIN);
					st = DRAINING_PRODUCERS;
				}
				break;

			case DRAINING_PRODUCERS:                             /* step 7 */
				if (all_in_state(PRODUCER, QUIESCENT))
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
	 * The round's statistics, once every worker is QUIESCENT: merge the producers' HyperLogLogs, reduce
	 * them (MAX) and the tallies (SUM) to rank 0, then zero everything round-scoped.  The reductions
	 * block; every rank reaches them after its drain (PROTOCOL.md §4.6).
	 */
	void end_round()
	{
		u64 sum[N_COUNTERS];
		snapshot(sum);
		u8 hll[HLL_REGISTERS] = {};
		for (int t = 0; t < params.n_threads; t++)
			if (ctx[t]->role == PRODUCER)
				hll_merge(hll, ctx[t]->hll.data());

		if (params.rank == 0) {
			MPI_Reduce(MPI_IN_PLACE, sum, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(MPI_IN_PLACE, hll, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
			controller.end_round(sum, hll);
		} else {
			MPI_Reduce(sum, NULL, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
			MPI_Reduce(hll, NULL, HLL_REGISTERS, MPI_UINT8_T, MPI_MAX, 0, params.mpi_comm);
		}

		for (int t = 0; t < params.n_threads; t++) {
			for (int k = 0; k < N_COUNTERS; k++)
				ctx[t]->ctr[k] = 0;
			for (int k = 0; k < (int) ctx[t]->hll.size(); k++)
				ctx[t]->hll[k] = 0;
		}
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
 * The one entry point of the engine.  Returns (i, x0, x1) --- the mixing function index and the two
 * colliding points --- or nothing if the search gave up after opts.max_versions rounds.  Threads pin
 * before they allocate: first touch places each object on its owner's NUMA node (PROTOCOL.md §4.1).
 */
template <class ProblemWrapper>
optional<tuple<u64,u64,u64>> run(const ProblemWrapper &wrapper, u64 nbytes_memory,
                                 const Options &opts, PRNG &prng)
{
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

		int want = params.place.thread_cpu[tid];
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
		int role = PRODUCER;
		if (tid == 0)
			role = COMM;
		else if (tid <= params.dicts_per_node)
			role = DICT;
		shared.ctx[tid] = std::make_unique<ThreadContext>(role, cpu, numa_node,
		                                                 params.place.thread_group[tid], params);
		if (role == DICT) {
			shared.shards[tid - 1] = std::make_unique<PcsDict>(params.jbits, params.w_shard);
			shared.coll_q[tid - 1] = std::make_unique<CollisionQueue>(params.coll_queue_capacity);
		}
		ThreadContext &me = *shared.ctx[tid];

		#pragma omp barrier             /* now we can inspect the shared context */

		if (tid == 0 && params.verbose) {
			std::vector<int> numa(params.n_threads), role(params.n_threads);
			for (int t = 0; t < params.n_threads; t++) {
				numa[t] = shared.ctx[t]->numa_node;
				role[t] = shared.ctx[t]->role;
			}
			params.place.report_measured(numa, role);
		}

		for (;;) {
			me.state = RUNNING;

			#pragma omp master
			comm.begin_round(wrapper.out_mask);

			#pragma omp barrier         /* makes the round header visible to all threads */

			if (shared.stop)
				break;

			if (me.role == COMM)
				comm.comm_round();          /* ends with end_round(), the statistics */
			else if (me.role == DICT)
				inserter_thread(me, params, shared, tid - 1);
			else
				walker_thread(me, wrapper, params, shared, tid - 1 - params.dicts_per_node);

			#pragma omp barrier         /* every thread is back.  Nothing depends on it */

			if (me.role == DICT)            /* later: a collective flush */
				shared.shards[tid - 1]->flush();
		}
	}

	optional<tuple<u64,u64,u64>> answer = comm.finish();

	/* detach waits for our last buffered sends; then restore the caller's buffer, if there was one */
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
