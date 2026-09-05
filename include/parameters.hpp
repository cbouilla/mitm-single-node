#ifndef MITM_PARAMETERS
#define MITM_PARAMETERS

#include <mpi.h>
#include <err.h>
#include <sched.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>

#include "tools.hpp"
#include "placement.hpp"

namespace mitm {

enum tags {TAG_POINTS, TAG_END_ROUND, TAG_REPORT, TAG_SOLUTION};

static constexpr int DP_WORDS = 2;     /* how many u64 per distinguished point on the wire */

struct DP {        /* one distinguished point (PROTOCOL.md §2.1) */
	u64 x;         /* the FULL endpoint */
	u64 jl;        /* the chain index j in the low jbits, the trail length above it */
};


/********************************** options **********************************/

struct Options {
	int inserters_per_node = 1;               /* dictionary shards per node. Users SHOULD set this themselves */
	int cache_level = 0;                      /* cache level a thread group sits in. 0 == auto-detect */
	double theta = -1;                        /* proportion of distinguished points. -1 == auto */
	bool verbose = true;                      /* print progress information (from rank 0 only) */

	/* --- below this line, the defaults should be just fine */

	/* nodes */
	MPI_Comm mpi_comm = MPI_COMM_WORLD;       /* one rank per node */

	/* thread layout inside a node */
	int walkers_per_node = 0;              /* 0 == fill the inherited affinity mask */
	bool bind_threads = true;              /* pin each thread to a CPU of the mask */

	/* algorithm */
	double alpha = 2.5;                       /* auto-chosen theta == alpha * sqrt(w/n) */
	double beta = 8;                          /* use each function variant for beta*w DPs */
	int dp_lenbits = 0;                       /* bits of trail length shipped with a DP.  0 == all that fit */
	u64 multiplier = 0x2545f4914f6cdd1dull;   /* to generate starting points */
	u64 max_versions = 0xffffffffffffffffull; /* how many functions to try before giving up */

	/* SPSC queue capacities, in DPs (rounded up to a power of two internally) */
	size_t walker_queue_capacity = 1024;   /* walker -> comm */
	size_t inserter_queue_capacity = 4096; /* comm -> inserter */

	/* control channel */
	int bsend_slack = 8;                   /* MPI_Bsend slots beyond the 2*n_nodes bounded traffic
	                                          (end-of-round signals + sentinels): our reports, our solution */

	/* bulk DP buffering */
	size_t buffer_capacity = 1500;         /* DPs per message (double-buffered per node) */
	int n_in_buffers = 8;                  /* posted MPI_Irecv slots (ANY_SOURCE) */

	/* collision queues, one per inserter */
	size_t coll_queue_capacity = 8192;     /* candidates buffered by each inserter for its group */
	size_t coll_per_chunk = 0;             /* candidates a walker retires per chunk; 0 == until its
	                                          batch is not full and the queue is empty */

	/* pacing */
	double ping_delay = 0.1;               /* max seconds between reports to the controller */
	int reports_per_round = 4;             /* also report every points_per_version/this DPs,
	                                          so a short round cannot overshoot beta*w */
	size_t chunk_size = 64;                /* vmixf iterations between queue / phase checks */

};


/********************************* parameters ********************************/

/* the rank and size of a communicator, so that Parameters can have them in its member init list */
static int comm_rank(MPI_Comm comm)
{
	int r;
	MPI_Comm_rank(comm, &r);
	return r;
}

static int comm_size(MPI_Comm comm)
{
	int n;
	MPI_Comm_size(comm, &n);
	return n;
}


/*
 * Everything the engine needs, derived once from (Options, RAM budget, n, m) at the top of run() and
 * const from then on: topology, thread layout and placement, dictionary size, difficulty.  Data only;
 * the Options are copied in, so their resolved values live here.
 */
struct Parameters : Options {
	/* MPI topology.  One rank per node, so `rank` IS the node index. */
	int rank;
	int n_nodes;                           /* == the number of MPI ranks */

	/* where this rank's threads go, and how many of them: placement.hpp */
	Placement place;

	/* thread layout */
	int n_threads;                         /* 1 + walkers_per_node + inserters_per_node */
	int n_walkers;                         /* total walker threads, all nodes */
	int n_inserters;                       /* total dictionary shards, all nodes */

	/* dictionary */
	u64 nbytes_memory;                     /* RAM budget per node */
	u64 w;                                 /* # slots in the (whole, distributed) dict */
	u64 w_shard;                           /* # slots per inserter thread */
	int jbits;                             /* bits of chain index stored with each DP */

	/* the second word of a DP: the chain index in the low jbits, the trail length above it */
	int lenbits;                           /* width of the length field */
	u64 len_sat;                           /* its largest value == "at least this, exact length unknown" */

	/* difficulty */
	bool theta_auto;                       /* theta was chosen from alpha */
	double auto_theta;                     /* what alpha chooses, whether or not it was used */
	u64 threshold;                         /* x is a DP iff x <= threshold */
	u64 dp_max_it;                         /* how many iterations to find a DP */
	u64 points_per_version;                /* #DP per version of the function */

	/*
	 * n, m are the wrapped problem's.  The thread layout is `place`'s, built first because it resolves
	 * `walkers_per_node` when the user left it at 0 (fill the affinity mask); everything here is the
	 * arithmetic that follows from it.
	 */
	Parameters(const Options &o, u64 nbytes_memory, int n, int m)
		: Options(o), rank(comm_rank(o.mpi_comm)), n_nodes(comm_size(o.mpi_comm)),
		  place(rank, o.inserters_per_node, o.walkers_per_node, o.bind_threads, o.cache_level),
		  nbytes_memory(nbytes_memory)
	{
		verbose = verbose && (rank == 0);
		walkers_per_node = place.n_walkers;

		n_threads = 1 + walkers_per_node + inserters_per_node;
		n_inserters = n_nodes * inserters_per_node;
		n_walkers = n_nodes * walkers_per_node;

		if (nbytes_memory == 0)
			errx(1, "the RAM budget per node must be given (and nonzero)");
		/* 8-byte slots, and a whole number of them per shard */
		w = (nbytes_memory * n_nodes) / sizeof(u64);
		w = (w / n_inserters) * n_inserters;
		if (w == 0)
			errx(1, "RAM budget too small: %" PRIu64 " bytes/node cannot hold one slot per shard",
			     nbytes_memory);
		assert(w % n_inserters == 0);
		w_shard = w / n_inserters;
		jbits = std::log2(10 * w) + 8;
		if (jbits > 56)
			errx(1, "dictionary too large: a %d-bit chain index leaves a slot no room for a length",
			     jbits);
		lenbits = dp_lenbits ? dp_lenbits : 64 - jbits;
		if (lenbits < 1 || jbits + lenbits > 64)
			errx(1, "--dp-len-bits %d does not fit beside a %d-bit chain index", lenbits, jbits);
		len_sat = make_mask(lenbits);

		auto_theta = alpha * std::sqrt(std::ldexp((double) w, -n));
		theta_auto = (theta < 0);
		if (theta_auto)
			theta = std::min(auto_theta, 1.0);
		threshold = pow(2, m) * theta;
		dp_max_it = 20 / theta;
		points_per_version = beta * w;
	}
};

}
#endif
