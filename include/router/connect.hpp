#ifndef MITM_ROUTER_CONNECT
#define MITM_ROUTER_CONNECT

#include "common.hpp"

namespace mitm {

/*
 * A sender gets its lines and its cache of blocks, a receiver its inbox, the service thread neither (an unused
 * inbox has one slot).  The object registers itself with the node, then holds the team's last barrier: this is
 * Router_Init's last stage, and the constructor runs nowhere else.  Past the barrier the service may read every
 * object, and thread 0 prints the measured layout.
 */
inline Router_thread::Router_thread(int role_, int index_, int group_, int domain_, int cpu_, int numa_,
                                    Router_node &rn)
	: swc(role_ == ROUTER_SENDER ? (Point *) router_alloc((size_t) rn.F * rn.swc_linesize * sizeof(Point)) : NULL),
	  swc_linesize(rn.swc_linesize), role(role_), index(index_),
	  global_id(role_ == ROUTER_SENDER ? rn.rank * rn.S + index_
	            : role_ == ROUTER_RECEIVER ? rn.rank * rn.R + index_ : -1),
	  group(group_), domain(domain_), cpu(cpu_), numa_node(numa_), node(rn), closed(0),
	  inbox(role_ == ROUTER_RECEIVER ? (size_t) rn.opt.inbox_blocks : 0)
{
	if (role == ROUTER_SENDER) {
		rn.senders[index] = this;
		rn.refill(*this);
	} else if (role == ROUTER_RECEIVER)
		rn.receivers[index] = this;
	#pragma omp barrier
	if (role == ROUTER_SERVICE && rn.rank == 0 && rn.opt.verbose)
		rn.plan.report_measured();
}

/* Thread 0's object owns the node.  The team's objects die in any order at the end of the region: a worker's
 * never touches the node here, and the node never touches a worker's. */
inline Router_thread::~Router_thread()
{
	::free(swc);
	if (role == ROUTER_SERVICE)
		delete &node;
}

/*
 * The shell, no MPI: the caller's communicator as is, `tag` for every message of ours -- the caller keeps both
 * and sends nothing else with that tag --, `lossy`, the options (NULL: the defaults) and room for the team's
 * roles.  Reached from Router_Init's single, on any thread of the team; connect() does the rest on thread 0.
 */
inline Router_node::Router_node(int n_threads, MPI_Comm mpi_comm, int tag_, bool lossy_, const Router_Opts *options)
	: comm(mpi_comm), tag(tag_), rank(-1), n_nodes(0), lossy(lossy_), opt(options ? *options : Router_Opts()),
	  input_closed(0)
{
	roles.assign(n_threads, -1);
	colors.assign(n_threads, ROUTER_GROUP_AUTO);
	masks.assign(n_threads, cpu_set_t());
}

/* Reached from thread 0's ~Router_thread, at the end of the region: only when nothing of ours is in flight, i.e.
 * quiescent or never used.  The workers' objects may be gone by now: nothing of theirs is read. */
inline Router_node::~Router_node()
{
	bool idle = phase == ROUTER_OPEN && ctr[ROUTER_BLOCKS] == 0 && ctr[ROUTER_MSGS_SENT] == 0;
	if (connected && not quiescent && not idle)
		errx(1, "Router: destroyed while not quiescent");
	for (size_t k = 0; k < in_req.size(); k++)
		if (in_req[k] != MPI_REQUEST_NULL) {
			MPI_Cancel(&in_req[k]);
			MPI_Wait(&in_req[k], MPI_STATUS_IGNORE);
		}
	free(dest);
	free(pool);
	free(n_valid);
	free(partial);
	free(blk_link);
}


/******************************** connection ********************************/

/*
 * Reached from Router_Init, thread 0 only, once every thread has declared its role and before any object is
 * built: the MPI checks, the agreement with the peers, the pool and the shared area, one block per destination,
 * the receives posted, a stash kept, every other block onto the free stack as one chain.
 */
inline void Router_node::connect()
{
	int provided;
	MPI_Query_thread(&provided);
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "Router: MPI_THREAD_FUNNELED is required (got %d)", provided);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &n_nodes);
	S = 0;
	R = 0;
	for (size_t t = 1; t < roles.size(); t++) {
		if (roles[t] == ROUTER_SENDER)
			S += 1;
		else
			R += 1;
	}
	if (S < 1 || R < 1)
		errx(1, "Router: a node needs at least one sender and one receiver (S=%d, R=%d)", S, R);

	/* every node must agree: lossy, S or R differing is fatal, an option differing means the defaults */
	const int N_AGREE = 13;
	double mine[N_AGREE] = {(double) lossy, (double) S, (double) R, (double) opt.block_points,
							(double) opt.swc_linesize, (double) opt.n_recv, (double) opt.inbox_blocks,
							(double) opt.sweep_blocks, (double) opt.credit, (double) opt.dests_per_node,
							(double) opt.pin, (double) opt.cache_level, (double) opt.group_size};
	double lo[N_AGREE], hi[N_AGREE];
	MPI_Allreduce(mine, lo, N_AGREE, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(mine, hi, N_AGREE, MPI_DOUBLE, MPI_MAX, comm);
	if (lo[0] != hi[0] || lo[1] != hi[1] || lo[2] != hi[2]) {
		if (rank == 0)
			warnx("Router: the nodes disagree on `lossy`, or have different numbers of senders or receivers");
		MPI_Abort(comm, 1);
	}
	for (int k = 3; k < N_AGREE; k++)
		if (lo[k] != hi[k]) {
			if (rank == 0)
				warnx("Router: the nodes disagree on the options; using the defaults everywhere");
			bool verbose = opt.verbose;
			opt = Router_Opts();
			opt.verbose = verbose;
			break;
		}

	if (not router_pow2(opt.block_points) || opt.block_points < 4)
		errx(1, "Router: block_points must be a power of two >= 4");
	if (opt.swc_linesize != 0 && (not router_pow2(opt.swc_linesize) || opt.swc_linesize < 4
								  || opt.swc_linesize > opt.block_points))
		errx(1, "Router: swc_linesize must be 0 or a power of two, 4 <= swc_linesize <= block_points");
	if (opt.swc_linesize != 0 && opt.block_points / opt.swc_linesize > ROUTER_MAX_L)
		errx(1, "Router: block_points / swc_linesize must be <= %u", ROUTER_MAX_L);
	if (opt.n_recv < 1 || opt.inbox_blocks < 1 || opt.sweep_blocks < 1 || opt.credit < 1)
		errx(1, "Router: n_recv, inbox_blocks, sweep_blocks and credit must be >= 1");
	if (opt.dests_per_node < 0)
		errx(1, "Router: dests_per_node must be >= 0");
	if (opt.group_size < 1)
		errx(1, "Router: group_size must be >= 1");
	if (opt.cache_level < 0)
		errx(1, "Router: cache_level must be >= 0");

	/* the union of the threads' own affinity masks: under OMP_PROC_BIND thread 0's alone would be one core */
	int n_threads = (int) roles.size();
	cpu_set_t umask;
	CPU_ZERO(&umask);
	for (int t = 0; t < n_threads; t++)
		for (int c = 0; c < CPU_SETSIZE; c++)
			if (CPU_ISSET(c, &masks[t]))
				CPU_SET(c, &umask);

	/* pinning must own the CPUs: refuse a mask too small for the team, or one shared with a co-hosted rank */
	if (opt.pin) {
		if (CPU_COUNT(&umask) < n_threads)
			errx(1, "Router: rank %d has %d CPUs in its affinity mask for %d threads; use --bind-to none,"
			        " --map-by slot:PE=n --bind-to core, or pin = false", rank, CPU_COUNT(&umask), n_threads);
		MPI_Comm host;
		MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &host);
		int hn = 0;
		int hr = 0;
		MPI_Comm_size(host, &hn);
		MPI_Comm_rank(host, &hr);
		if (hn > 1) {
			std::vector<cpu_set_t> all(hn);
			MPI_Allgather(&umask, sizeof(cpu_set_t), MPI_BYTE, all.data(), sizeof(cpu_set_t), MPI_BYTE, host);
			for (int o = 0; o < hn; o++) {
				if (o == hr)
					continue;
				int inter = 0;
				for (int c = 0; c < CPU_SETSIZE; c++)
					if (CPU_ISSET(c, &umask) && CPU_ISSET(c, &all[o]))
						inter += 1;
				if (inter > 0) {
					warnx("Router: rank %d shares CPUs with another rank on its host; use --map-by numa"
					      " --bind-to numa, --map-by slot:PE=n, or pin = false", rank);
					MPI_Abort(comm, 1);
				}
			}
		}
		MPI_Comm_free(&host);
	}

	plan = RouterPlacement(umask, roles, colors, opt.pin, opt.cache_level, opt.group_size, rank);
	if (rank == 0 && opt.verbose)
		plan.report();
	if (opt.pin && plan.thread_cpu[0] >= 0 && pin_to_cpu(plan.thread_cpu[0]) < 0)
		warn("Router: rank %d: cannot pin the service thread to CPU %d", rank, plan.thread_cpu[0]);

	int P = n_nodes;
	credit_used.assign(P, 0);
	seq_sent.assign(P, 0);
	end_msg.assign(P, RouterMsgHdr());
	end_req.assign(P, MPI_REQUEST_NULL);
	end_sent.assign(P, 0);
	n_data_recv.assign(P, 0);
	end_seq.assign(P, -1);

	per_node = opt.dests_per_node ? opt.dests_per_node : R;
	if (per_node < R)
		errx(1, "Router: dests_per_node (%d) must be at least the receivers per node (%d)", per_node, R);
	F = per_node * n_nodes;

	size_t swc = opt.swc_linesize;
	if (swc == 0) {                       /* the lines fit L2, and a block switch is at most one line in 8 */
		size_t budget = (512 << 10) / ((size_t) F * sizeof(Point));
		swc = 4;
		while (swc * 2 <= budget && swc * 2 <= opt.block_points / 8)
			swc *= 2;
	}
	while (opt.block_points / swc > ROUTER_MAX_L)
		swc *= 2;
	swc_linesize = swc;
	L = (u32) (opt.block_points / swc);
	nv_stride = (L + 63) & ~(u32) 63;
	block_bytes = ROUTER_HDR_BYTES + opt.block_points * sizeof(Point);

	/* installed, the stash, the slots, what the receivers may hold, the senders' caches, pending, and sealed or
	 * parked */
	size_t slack = (size_t) F > 32 * (size_t) S ? (size_t) F : 32 * (size_t) S;
	size_t nb = (size_t) F + 2 * ROUTER_BATCH + 2 * (size_t) opt.n_recv
				+ (size_t) R * ((size_t) opt.inbox_blocks + 1) + (size_t) S * (ROUTER_BATCH + 1) + slack;
	if (nb >= ROUTER_NONE)
		errx(1, "Router: too many blocks (%zu)", nb);
	n_blocks = (u32) nb;
	dest = (std::atomic<u64> *) router_alloc((size_t) F * ROUTER_DEST_WORDS * sizeof(u64));
	pool = (char *) router_alloc((size_t) n_blocks * block_bytes);
	n_valid = (std::atomic<u8> *) router_alloc((size_t) n_blocks * nv_stride);
	partial_cap = ((size_t) S * (swc_linesize - 1) + 3) & ~(size_t) 3;
	partial = (Point *) router_alloc((size_t) F * partial_cap * sizeof(Point));
	blk_link = (std::atomic<u64> *) router_alloc((size_t) n_blocks * sizeof(u64));
	free_list.reserve(n_blocks);
	for (u32 b = 0; b < n_blocks; b++) {
		blk_link[b].store(ROUTER_NONE, std::memory_order_relaxed);
		free_list.push_back(b);
	}
	blk_dest.assign(n_blocks, -1);
	blk_count.assign(n_blocks, 0);
	park_head.assign(R + n_nodes, ROUTER_NONE);
	park_tail.assign(R + n_nodes, ROUTER_NONE);
	pending.reserve(S + 1);
	inbox_pushed.assign(R, 0);

	out_req.assign(opt.n_recv, MPI_REQUEST_NULL);
	out_blk.assign(opt.n_recv, ROUTER_NONE);
	out_peer.assign(opt.n_recv, -1);
	out_free.clear();
	for (int k = 0; k < opt.n_recv; k++)
		out_free.push_back(k);
	in_req.assign(opt.n_recv, MPI_REQUEST_NULL);
	in_blk.assign(opt.n_recv, ROUTER_NONE);
	mpi_idx.assign(opt.n_recv, 0);
	mpi_st.assign(opt.n_recv, MPI_Status());
	senders.assign(S, NULL);
	receivers.assign(R, NULL);

	for (int d = 0; d < F; d++)
		dest[ROUTER_DEST_WORDS * d].store((u64) stash_pop(), std::memory_order_relaxed);
	for (int k = 0; k < opt.n_recv; k++)
		repost(k);
	if (free_list.size() > ROUTER_BATCH)
		spill(free_list.size() - ROUTER_BATCH);
	connected = true;
	if (rank == 0 && opt.verbose)
		banner();
}

/* the sizes, once.  Reached from connect(), on rank 0 only, when opt.verbose. */
inline void Router_node::banner() const
{
	double blocks = (double) n_blocks * block_bytes;
	double lines = (double) S * F * swc_linesize * sizeof(Point);
	double inboxes = (double) R * opt.inbox_blocks * sizeof(Point);
	double closing = (double) F * partial_cap * sizeof(Point);
	printf("Router: %s, %d node%s, %d sender%s and %d receiver%s per node, %d destinations (%d per node)\n",
	       lossy ? "lossy" : "lossless", n_nodes, n_nodes > 1 ? "s" : "", S, S > 1 ? "s" : "",
	       R, R > 1 ? "s" : "", F, per_node);
	printf("Router: blocks of %zu points, lines of %zu (%u lines per block), %u blocks (%.1f MB, %u cached per"
	       " sender), messages of %zu bytes\n", opt.block_points, swc_linesize, L, n_blocks, blocks / 1e6,
	       ROUTER_BATCH, block_bytes);
	printf("Router: private lines %.1f MB (%.1f MB per sender), closing buffers %.1f MB (%zu points per destination)\n",
	       lines / 1e6, lines / 1e6 / S, closing / 1e6, partial_cap);
	printf("Router: inboxes of %d blocks (%.1f MB), %d receive and %d send slots, credit %d per peer, stash batches"
	       " of %u, %d sealed blocks per turn\n", opt.inbox_blocks, inboxes / 1e6, opt.n_recv, opt.n_recv, opt.credit,
	       ROUTER_BATCH, opt.sweep_blocks);
}

/*
 * Collective over every thread of the team and every node, the same arguments on every thread except `group`,
 * a thread's own (ROUTER_GROUP_AUTO to let the Router form the groups, or a color in 0..G-1).  Thread 0 is the
 * service whatever `role` says.  Each thread records its role, its group and its own affinity mask; thread 0
 * connects and builds the plan; then every thread pins itself, confirms the kernel obeyed, and builds its
 * object, which every later call of this thread takes; the caller keeps it, named, for the team's whole life.
 */
[[nodiscard]] inline Router_thread Router_Init(int role, int group, MPI_Comm comm, int tag, bool lossy,
                                               const Router_Opts *opts)
{
	int tid = omp_get_thread_num();
	Router_node *rn = NULL;
	#pragma omp single copyprivate(rn)
	rn = new Router_node(omp_get_num_threads(), comm, tag, lossy, opts);
	if (tid == 0)
		rn->roles[0] = ROUTER_SERVICE;
	else if (role == ROUTER_SENDER || role == ROUTER_RECEIVER)
		rn->roles[tid] = role;
	else
		errx(1, "Router_Init: role must be ROUTER_SENDER or ROUTER_RECEIVER");
	rn->colors[tid] = (tid == 0) ? ROUTER_GROUP_AUTO : group;
	if (sched_getaffinity(0, sizeof(cpu_set_t), &rn->masks[tid]) != 0)
		err(1, "Router_Init: sched_getaffinity");
	#pragma omp barrier
	if (tid == 0)
		rn->connect();
	#pragma omp barrier

	int want = rn->plan.thread_cpu[tid];  /* pin as the plan says, then confirm the kernel put us there */
	if (want >= 0 && pin_to_cpu(want) < 0) {
		warn("Router: rank %d: cannot pin thread %d to CPU %d", rn->rank, tid, want);
		#pragma omp atomic
		rn->misplaced += 1;
	}
	unsigned cpu = 0;
	unsigned numa = 0;
	if (syscall(SYS_getcpu, &cpu, &numa, NULL) != 0)
		err(1, "Router_Init: getcpu");
	if (want >= 0 && (int) cpu != want) {
		warnx("Router: rank %d: thread %d asked for CPU %d, runs on CPU %u", rn->rank, tid, want, cpu);
		#pragma omp atomic
		rn->misplaced += 1;
	}
	#pragma omp barrier                   /* everybody is pinned, or nobody builds an object */
	if (tid == 0 && rn->misplaced > 0)
		MPI_Abort(comm, 1);

	int mine = rn->roles[tid];
	int index = 0;                        /* among the threads of my role: those before me */
	for (int t = 0; t < tid; t++)
		if (rn->roles[t] == mine)
			index += 1;
	return Router_thread(mine, index, rn->plan.thread_group[tid], rn->plan.thread_domain[tid],
	                     (int) cpu, (int) numa, *rn);   /* registers itself, then holds the team's last barrier */
}

/* Any thread.  What the team looks like, once connected: the senders and receivers of this node and of all. */
inline int Router_local_num_send(const Router_thread &rt)
{
	return rt.node.S;
}

inline int Router_local_num_recv(const Router_thread &rt)
{
	return rt.node.R;
}

inline int Router_num_send(const Router_thread &rt)
{
	return rt.node.S * rt.node.n_nodes;
}

inline int Router_num_recv(const Router_thread &rt)
{
	return rt.node.R * rt.node.n_nodes;
}

/* Sender or receiver.  The caller's rank among the threads of its role on this node, from 0. */
inline int Router_local_rank(const Router_thread &rt)
{
	if (rt.role != ROUTER_SENDER && rt.role != ROUTER_RECEIVER)
		errx(1, "Router_local_rank: not a sender or receiver thread");
	return rt.index;
}

/* Sender or receiver.  The caller's rank among the threads of its role on all nodes, from 0: the nodes' in
 * order, each node's by local rank.  A receiver's is the `dest` that reaches it. */
inline int Router_rank(const Router_thread &rt)
{
	if (rt.role != ROUTER_SENDER && rt.role != ROUTER_RECEIVER)
		errx(1, "Router_rank: not a sender or receiver thread");
	return rt.global_id;
}

/* Sender or receiver.  Its group, 0..Router_num_groups(rt)-1; -1 for the service thread. */
inline int Router_group(const Router_thread &rt)
{
	return rt.group;
}

/* Any thread.  The groups this node was cut into. */
inline int Router_num_groups(const Router_thread &rt)
{
	return rt.node.plan.n_groups;
}

/* Sender or receiver.  The receivers in the caller's group: a producer pairs with one of them. */
inline int Router_group_num_receivers(const Router_thread &rt)
{
	if (rt.group < 0)
		errx(1, "Router_group_num_receivers: not a sender or receiver thread");
	return (int) rt.node.plan.group_receivers[rt.group].size();
}

/* Sender or receiver.  The local index (0..R-1) of the i-th receiver in the caller's group. */
inline int Router_group_receiver(int i, const Router_thread &rt)
{
	if (rt.group < 0)
		errx(1, "Router_group_receiver: not a sender or receiver thread");
	return rt.node.plan.group_receivers[rt.group].at(i);
}

/* Any thread.  The cache domain it was pinned in; -1 when the Router did not pin. */
inline int Router_domain(const Router_thread &rt)
{
	return rt.domain;
}

/* Any thread.  This node's cache domains holding a CPU of the mask; 0 when not pinned. */
inline int Router_num_domains(const Router_thread &rt)
{
	return rt.node.plan.n_domains;
}

/* Any thread.  The CPU and NUMA node the kernel reports for it (momentary when not pinned). */
inline int Router_cpu(const Router_thread &rt)
{
	return rt.cpu;
}

inline int Router_numa_node(const Router_thread &rt)
{
	return rt.numa_node;
}

}
#endif
