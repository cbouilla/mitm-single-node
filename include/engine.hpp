#ifndef MITM_ENGINE
#define MITM_ENGINE

#include <cassert>
#include <cstring>
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
#include "controller.hpp"

namespace mitm {


/******************************* the comm thread ******************************/

/*
 * Thread 0's object: all of the rank's MPI.  Routing, the receive pool, the control channel, the
 * round's state machine (comm_round), the end-of-round reductions and, on rank 0, the controller.
 * Built by run() on the main thread, before the team.  The scheme supplies its parameters, its
 * counters and its overflow policy; the routing is the same for every scheme.  PROTOCOL.md §2 and §4.
 */
template <class Scheme>
class CommThread {
public:
	using Params = typename Scheme::Params;
	static constexpr int N = Scheme::N_COUNTERS;

	const Params &params;
	PRNG &prng;
	SharedContext<Scheme> &shared;
	std::vector<std::unique_ptr<ThreadContext<Scheme>>> &ctx;   /* == shared.ctx: n_threads, indexed by tid */

	OutBuffers outbuf;

	/* the receive pool: n_in_buffers always-posted ANY_SOURCE receives, sized by that knob and not by the peers */
	size_t in_cap;                                        /* u64 per buffer */
	std::vector<std::vector<u64>> in_data;
	std::vector<MPI_Request> in_req;
	std::vector<int> in_done_idx;                         /* MPI_Testsome output */
	std::vector<MPI_Status> in_done_st;                   /* MPI_Testsome output */
	int n_sentinels = 0;                                  /* nodes done sending this round */

	MPI_Request req_end_round = MPI_REQUEST_NULL;         /* the TAG_END_ROUND receive, always posted */
	Controller<Scheme> controller;

	bool golden_sent = false;
	bool round_over = false;           /* our TAG_END_ROUND has arrived */

	/*
	 * A lossless scheme's backlog (Scheme::LOSSLESS, PROTOCOL.md §5): what found its channel full waits
	 * here and is retried first on the next turn, instead of being dropped.  The carry holds at most one
	 * batch, because nothing new is popped while it is not empty; a received buffer waits in `in_ready`
	 * until every point of it is delivered, and is reposted only then.  A dropping scheme leaves all of
	 * this empty, so the drain tests that read it (§4.5, steps 2 and 5) cost it nothing.
	 */
	std::vector<Point> carry;          /* points popped from a producer ring that found their channel full */
	size_t carry_n = 0;
	std::vector<int> in_ready;         /* received buffers not yet fully delivered: a ring of indices into in_data */
	std::vector<int> in_count;         /* words received into each buffer */
	size_t in_ready_head = 0;          /* the oldest */
	size_t in_ready_n = 0;
	size_t in_off = 0;                 /* how far the oldest has been delivered, in words */

	CommThread(const Params &params, PRNG &prng, SharedContext<Scheme> &shared)
		: params(params), prng(prng), shared(shared), ctx(shared.ctx),
		  outbuf(params.mpi_comm, params.n_nodes, params.buffer_capacity),
		  in_cap(POINT_WORDS * params.buffer_capacity), in_data(params.n_in_buffers),
		  in_req(params.n_in_buffers, MPI_REQUEST_NULL),
		  in_done_idx(params.n_in_buffers), in_done_st(params.n_in_buffers),
		  controller(params), carry(64), in_ready(params.n_in_buffers), in_count(params.n_in_buffers, 0)
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
	 * ping-pong between the two cores once per point (PROBLEM.md §3, §5.C).  For a dropping scheme a
	 * bucket never outlives the call that filled it -- route_producer_queues() and poll_incoming() both
	 * flush before returning; a lossless one keeps what did not fit in the ring, and the drain waits for
	 * the buckets to empty (PROTOCOL.md §3.1, §4.5).
	 */
	static constexpr size_t BUCKET_CAP = 128;             /* points staged per shard */
	std::vector<Point> bucket;                            /* dicts_per_node runs of BUCKET_CAP */
	std::vector<size_t> bucket_n;                         /* points staged for each shard */

	/* hand shard `slot`'s staged run to its dict thread.  What did not fit is dropped, or kept for the next turn */
	void flush_bucket(int slot)
	{
		size_t n = bucket_n[slot];
		if (n == 0)
			return;
		Point *run = &bucket[slot * BUCKET_CAP];
		size_t k = ctx[1 + slot]->q->push_bulk(run, n);          /* dict thread `slot` is thread 1 + slot */
		if constexpr (Scheme::LOSSLESS) {
			for (size_t t = k; t < n; t++)
				run[t - k] = run[t];
			bucket_n[slot] = n - k;
		} else {
			ctx[0]->ctr[Scheme::DROP_DICTQ] += n - k;
			bucket_n[slot] = 0;
		}
	}

	void flush_local()
	{
		for (int slot = 0; slot < params.dicts_per_node; slot++)
			flush_bucket(slot);
	}

	/* stage a point for the dict thread that owns its shard, on this node.  false == bucket and ring both full */
	bool deliver_local(const Point &p)
	{
		int slot = (int) ((p.key / params.n_nodes) % params.dicts_per_node);
		if (bucket_n[slot] == BUCKET_CAP) {
			flush_bucket(slot);
			if (bucket_n[slot] == BUCKET_CAP)
				return false;
		}
		bucket[slot * BUCKET_CAP + bucket_n[slot]] = p;
		bucket_n[slot] += 1;
		return true;
	}

	/* route one point, to a local bucket or to an outgoing buffer.  false == its channel is full */
	bool place(const Point &p)
	{
		int dst = (int) (p.key % params.n_nodes);
		if (dst == params.rank)
			return deliver_local(p);
		return outbuf.push(p, dst);
	}

	/* is anything received or staged still waiting for a dict ring?  Always false for a dropping scheme */
	bool local_pending()
	{
		if (in_ready_n > 0)
			return true;
		for (int slot = 0; slot < params.dicts_per_node; slot++)
			if (bucket_n[slot] > 0)
				return true;
		return false;
	}

	/*
	 * Pull whatever the producer threads have produced and route it.  A lossless scheme first retries
	 * what the last turn could not place -- the held buckets, then the carry -- and pops nothing new
	 * while the carry is not empty; a batch that stalls ends the sweep and the other rings wait a turn.
	 */
	bool route_producer_queues()
	{
		if constexpr (Scheme::LOSSLESS) {
			flush_local();
			size_t kept = 0;
			for (size_t t = 0; t < carry_n; t++)
				if (not place(carry[t]))
					carry[kept++] = carry[t];
			carry_n = kept;
			if (carry_n > 0) {
				ctx[0]->ctr[Scheme::STALL_OUT] += 1;
				return false;
			}
		}
		Point batch[64];
		bool moved = false;
		for (int s = 0; s < params.producers_per_node; s++) {
			int t = 1 + params.dicts_per_node + s;              /* producer s is thread t */
			size_t k = ctx[t]->q->pop_bulk(batch, 64);
			if (k > 0)
				moved = true;
			for (size_t u = 0; u < k; u++) {
				if (place(batch[u]))
					continue;
				if constexpr (Scheme::LOSSLESS)
					carry[carry_n++] = batch[u];
				else
					ctx[0]->ctr[Scheme::DROP_OUT] += 1;   /* a bucket always takes a point here: a remote refusal */
			}
			if constexpr (Scheme::LOSSLESS)
				if (carry_n > 0)
					break;
		}
		flush_local();
		return moved;
	}

	/* post buffer k's receive again */
	void repost(int k)
	{
		MPI_Irecv(in_data[k].data(), in_cap, MPI_UINT64_T, MPI_ANY_SOURCE, TAG_POINTS, params.mpi_comm, &in_req[k]);
	}

	/*
	 * Take delivery of whatever arrived over MPI, oldest buffer first, and repost each receive once its
	 * buffer is fully delivered.  A sentinel is counted at once.  A lossless scheme's delivery stops at
	 * a point whose dict ring is full and resumes there on the next turn, which is why the WAITING exit
	 * also asks for every received buffer to be delivered (PROTOCOL.md §4.5, §5).
	 */
	void poll_incoming()
	{
		int outcount = 0;
		MPI_Testsome(in_req.size(), in_req.data(), &outcount, in_done_idx.data(), in_done_st.data());
		if (outcount != MPI_UNDEFINED)
			for (int t = 0; t < outcount; t++) {
				int k = in_done_idx[t];
				MPI_Get_count(&in_done_st[t], MPI_UINT64_T, &in_count[k]);
				if (in_count[k] == 0) {
					n_sentinels += 1;             /* that node has finished the round */
					repost(k);
				} else {
					in_ready[(in_ready_head + in_ready_n) % in_ready.size()] = k;
					in_ready_n += 1;
				}
			}

		while (in_ready_n > 0) {
			int k = in_ready[in_ready_head];
			bool stalled = false;
			for (; in_off + POINT_WORDS <= (size_t) in_count[k]; in_off += POINT_WORDS) {
				Point p = {in_data[k][in_off], in_data[k][in_off + 1]};
				if (not deliver_local(p)) {
					stalled = true;
					break;
				}
			}
			if (stalled)
				break;
			repost(k);
			in_ready_head = (in_ready_head + 1) % in_ready.size();
			in_ready_n -= 1;
			in_off = 0;
		}
		flush_local();
		if constexpr (Scheme::LOSSLESS)
			if (local_pending())
				ctx[0]->ctr[Scheme::STALL_IN] += 1;
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
	void snapshot(u64 sum[N])
	{
		#pragma omp flush
		for (int k = 0; k < N; k++)
			sum[k] = 0;
		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N; k++)
				sum[k] += ctx[t]->ctr[k];
	}

	/*
	 * The round header: the scheme's words plus `stop`, drawn by rank 0 from the previous header and
	 * the controller's verdict, broadcast, stored in the SharedContext (PROTOCOL.md §4.2).
	 */
	template <class Wrapper>
	void begin_round(const Wrapper &wrapper)
	{
		constexpr int HW = sizeof(typename Scheme::Header) / sizeof(u64);
		u64 msg[HW + 1];
		if (params.rank == 0) {
			u64 stop = controller.stop;
			Scheme::next_header(params, wrapper, prng, shared.header, stop);
			memcpy(msg, &shared.header, sizeof(shared.header));
			msg[HW] = stop;
		}
		MPI_Bcast(msg, HW + 1, MPI_UINT64_T, 0, params.mpi_comm);
		memcpy(&shared.header, msg, sizeof(shared.header));
		shared.stop = msg[HW];

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

	/* step 1 of the drain: the producers still RUNNING are told to HOLD; one already HELD (exhausted) is left be */
	void hold_producers()
	{
		for (int t = 0; t < params.n_threads; t++)
			if (ctx[t]->role == PRODUCER && ctx[t]->state.load(std::memory_order_acquire) == RUNNING)
				ctx[t]->state.store(HOLD, std::memory_order_release);
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

		u64 prev[N] = {0};                 /* the tallies at the previous report */
		u64 poll_tick = 0;

		for (;;) {
			bool moved = route_producer_queues();
			poll_incoming();
			service_control();

			if (not golden_sent && shared.golden_found.load(std::memory_order_acquire)) {
				MPI_Bsend(shared.golden, SOL_NWORDS, MPI_UINT64_T, 0, TAG_SOLUTION, params.mpi_comm);
				golden_sent = true;
			}

			/* report by time, and by volume too: on a timer alone a short round overshoots (§2.2) */
			if ((++poll_tick & 0xff) == 0) {
				u64 cur[N];
				snapshot(cur);
				bool due = cur[Scheme::PACING] - prev[Scheme::PACING] >= params.report_points
				        || wtime() - last_ping >= params.ping_delay;
				if (due) {
					last_ping = wtime();
					u64 msg[N];
					for (int k = 0; k < N; k++) {
						msg[k] = cur[k] - prev[k];
						prev[k] = cur[k];
					}
					MPI_Bsend(msg, N, MPI_UINT64_T, 0, TAG_REPORT, params.mpi_comm);
				}
			}

			/* our phase: nobody else reads it, hence relaxed */
			int st = ctx[0]->state.load(std::memory_order_relaxed);
			assert(st == RUNNING || st == COLLECTING || not moved);      /* held producers produce nothing */

			switch (st) {
			case RUNNING:                                        /* step 1: told to, or every producer exhausted */
				if (round_over || all_in_state(PRODUCER, HELD)) {
					hold_producers();
					st = COLLECTING;
				}
				break;

			case COLLECTING:                                     /* step 2 */
				if (all_in_state(PRODUCER, HELD) && not moved && producer_queues_empty() && carry_n == 0)
					st = FLUSHING;
				break;

			case FLUSHING:                                       /* steps 3, 4 */
				if (outbuf.flush_poll()) {
					outbuf.send_sentinels();
					st = WAITING;
				}
				break;

			case WAITING:                                        /* steps 5, 6 */
				if (n_sentinels == params.n_nodes && not local_pending()) {
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
	 * The round's statistics, once every worker is QUIESCENT: the tallies (SUM) and the scheme's own
	 * statistics to rank 0, then every node's golden slot, so that a pair found during the drain is
	 * never lost (PROTOCOL.md §4.6); then zero everything round-scoped.  The collectives block; every
	 * rank reaches them after its drain.
	 */
	void end_round()
	{
		u64 sum[N];
		snapshot(sum);
		typename Scheme::RoundStats stats;
		stats.collect(ctx);                /* the producers' own statistics, merged and zeroed */

		if (params.rank == 0)
			MPI_Reduce(MPI_IN_PLACE, sum, N, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
		else
			MPI_Reduce(sum, NULL, N, MPI_UINT64_T, MPI_SUM, 0, params.mpi_comm);
		stats.reduce(params.mpi_comm, params.rank);

		/* the answer, exactly: (found, i, x0, x1) from every node; rank 0 keeps the first one found */
		constexpr int GW = 1 + SOL_NWORDS;
		u64 mine[GW] = {shared.golden_found.load(std::memory_order_acquire),
		                shared.golden[SOL_I], shared.golden[SOL_X0], shared.golden[SOL_X1]};
		std::vector<u64> every((params.rank == 0) ? GW * params.n_nodes : 0);
		MPI_Gather(mine, GW, MPI_UINT64_T, every.data(), GW, MPI_UINT64_T, 0, params.mpi_comm);

		if (params.rank == 0) {
			for (int r = 0; r < params.n_nodes; r++)
				if (every[GW * r])
					controller.record(every[GW * r + 1 + SOL_I], every[GW * r + 1 + SOL_X0],
					                  every[GW * r + 1 + SOL_X1]);
			controller.end_round(sum, stats);
		}

		for (int t = 0; t < params.n_threads; t++)
			for (int k = 0; k < N; k++)
				ctx[t]->ctr[k] = 0;
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
 * The one entry point of the engine core, for any scheme.  Returns (i, x0, x1) --- the scheme's round
 * word and the two colliding points --- or nothing if the search stopped without one.  Threads pin
 * before they allocate: first touch places each object on its owner's NUMA node (PROTOCOL.md §4.1).
 */
template <class Scheme, class Wrapper>
optional<tuple<u64,u64,u64>> run(const Wrapper &wrapper, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI: MPI_THREAD_FUNNELED is required (got %d).  Use MPI_Init_thread.", provided);

	const typename Scheme::Params params(opts, nbytes_memory, wrapper.n, wrapper.m);

	/* The MPI_Bsend buffer of the engine (preserving the previous one) */
	static_assert((int) Scheme::N_COUNTERS >= (int) SOL_NWORDS, "the Bsend slots are sized for a report");
	void *prev_bsend_buf = NULL;
	int prev_bsend_size = 0;
	MPI_Buffer_detach(&prev_bsend_buf, &prev_bsend_size);
	size_t slots = 2 * (size_t) params.n_nodes + params.bsend_slack;
	size_t msg = Scheme::N_COUNTERS * sizeof(u64) + MPI_BSEND_OVERHEAD;
	std::vector<char> bsend_buf(slots * msg);
	MPI_Buffer_attach(bsend_buf.data(), (int) bsend_buf.size());

	/* safety check: all ranks evaluate the same function (no uninitialized state) */
	u64 test[3];
	test[0] = prng.rand();
	test[1] = prng.rand();
	test[2] = wrapper.self_test(test[0], test[1]);
	MPI_Bcast(test, 3, MPI_UINT64_T, 0, params.mpi_comm);
	assert(test[2] == wrapper.self_test(test[0], test[1]));

	if (params.verbose)
		Scheme::banner(params, prng.seed);

	SharedContext<Scheme> shared(params);
	CommThread<Scheme> comm(params, prng, shared);
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
		shared.ctx[tid] = std::make_unique<ThreadContext<Scheme>>(role, cpu, numa_node,
		                                                         params.place.thread_group[tid], params);
		if (role == DICT)
			Scheme::build_dict(shared, params, tid - 1);
		ThreadContext<Scheme> &me = *shared.ctx[tid];

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
			comm.begin_round(wrapper);

			#pragma omp barrier         /* makes the round header visible to all threads */

			if (shared.stop)
				break;
			/* this thread's copy of the round's header: thread 0 rewrites shared.header at the next round's
			   start, while a slower thread may still be in after_round() */
			const typename Scheme::Header header = shared.header;

			if (me.role == COMM)
				comm.comm_round();          /* ends with end_round(), the statistics */
			else if (me.role == DICT)
				Scheme::dict_thread(me, wrapper, params, shared, tid - 1);
			else
				Scheme::producer_thread(me, wrapper, params, shared, tid - 1 - params.dicts_per_node);

			#pragma omp barrier         /* every thread is back.  Nothing depends on it */

			if (me.role == DICT)
				Scheme::after_round(shared, params, tid - 1, header);
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
