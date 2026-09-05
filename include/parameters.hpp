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

static constexpr int POINT_WORDS = 2;     /* how many u64 per point on the wire */

/*
 * The unit of traffic, from a producer to the shard that owns it (PROTOCOL.md §2.1): the routing key
 * in full and one word of payload.  What the words mean is the scheme's business -- PCS ships a
 * distinguished point as (endpoint, chain index | trail length), the direct scheme (image, preimage).
 */
struct Point {
	u64 key;       /* routed on it: node = key % n_nodes, shard = (key / n_nodes) % dicts_per_node */
	u64 val;       /* the payload */
};


/********************************** options **********************************/

/* Every user knob of every scheme, all with a working default.  A scheme ignores the other's. */
struct Options {
	int dicts_per_node = 1;                   /* dictionary shards per node. Users SHOULD set this themselves */
	int cache_level = 0;                      /* cache level a thread group sits in. 0 == auto-detect */
	double theta = -1;                        /* PCS: proportion of distinguished points. -1 == auto */
	bool verbose = true;                      /* print progress information (from rank 0 only) */

	/* --- below this line, the defaults should be just fine */

	/* nodes */
	MPI_Comm mpi_comm = MPI_COMM_WORLD;       /* one rank per node */

	/* thread layout inside a node */
	int producers_per_node = 0;            /* 0 == fill the inherited affinity mask */
	bool bind_threads = true;              /* pin each thread to a CPU of the mask */

	/* PCS */
	double alpha = 2.5;                       /* auto-chosen theta == alpha * sqrt(w/n) */
	double beta = 8;                          /* use each function variant for beta*w DPs */
	int dp_lenbits = 0;                       /* bits of trail length shipped with a DP.  0 == all that fit */
	u64 multiplier = 0x2545f4914f6cdd1dull;   /* to generate starting points */
	u64 max_versions = 0xffffffffffffffffull; /* how many functions to try before giving up */

	/* SPSC queue capacities, in points (rounded up to a power of two internally) */
	size_t producer_queue_capacity = 1024; /* producer -> comm */
	size_t dict_queue_capacity = 4096;     /* comm -> dict thread */

	/* control channel */
	int bsend_slack = 8;                   /* MPI_Bsend slots beyond the 2*n_nodes bounded traffic
	                                          (end-of-round signals + sentinels): our reports, our solution */

	/* bulk point buffering */
	size_t buffer_capacity = 1500;         /* points per message (double-buffered per node) */
	int n_in_buffers = 8;                  /* posted MPI_Irecv slots (ANY_SOURCE) */

	/* PCS: collision queues, one per dict thread */
	size_t coll_queue_capacity = 8192;     /* candidates buffered by each dict thread for its group */
	size_t coll_per_chunk = 0;             /* candidates a producer retires per chunk; 0 == until its
	                                          batch is not full and the queue is empty */

	/* pacing */
	double ping_delay = 0.1;               /* max seconds between reports to the controller */
	int reports_per_round = 4;             /* also report every round_points/this points, so a short
	                                          round cannot overshoot its volume */
	size_t chunk_size = 64;                /* PCS: vmixf iterations between queue / phase checks */

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
 * What the engine core needs, derived once from (Options, RAM budget) at the top of run() and const
 * from then on: topology, thread layout and placement, the dictionary's slots.  Data only; the Options
 * are copied in, so their resolved values live here.  A scheme derives its own Params from it.
 */
struct Parameters : Options {
	/* MPI topology.  One rank per node, so `rank` IS the node index. */
	int rank;
	int n_nodes;                           /* == the number of MPI ranks */

	/* where this rank's threads go, and how many of them: placement.hpp */
	Placement place;

	/* thread layout */
	int n_threads;                         /* 1 + producers_per_node + dicts_per_node */
	int n_producers;                       /* total producer threads, all nodes */
	int n_dicts;                           /* total dictionary shards, all nodes */

	/* dictionary: 8-byte slots, a whole number of them per shard */
	u64 nbytes_memory;                     /* RAM budget per node */
	u64 w;                                 /* # slots in the (whole, distributed) dictionary */
	u64 w_shard;                           /* # slots per dict thread */

	/* report pacing: also report every this many points (PROTOCOL.md §2.2).  The scheme sets it */
	u64 report_points = 1;

	/*
	 * The thread layout is `place`'s, built first because it resolves `producers_per_node` when the
	 * user left it at 0 (fill the affinity mask); everything here is the arithmetic that follows.
	 * `producer_per_dict` is the scheme's: does every dict thread need a producer of its own group?
	 */
	Parameters(const Options &o, u64 nbytes_memory, bool producer_per_dict)
		: Options(o), rank(comm_rank(o.mpi_comm)), n_nodes(comm_size(o.mpi_comm)),
		  place(rank, o.dicts_per_node, o.producers_per_node, o.bind_threads, o.cache_level, producer_per_dict),
		  nbytes_memory(nbytes_memory)
	{
		verbose = verbose && (rank == 0);
		producers_per_node = place.n_producers;

		n_threads = 1 + producers_per_node + dicts_per_node;
		n_dicts = n_nodes * dicts_per_node;
		n_producers = n_nodes * producers_per_node;

		if (nbytes_memory == 0)
			errx(1, "the RAM budget per node must be given (and nonzero)");
		w = (nbytes_memory * n_nodes) / sizeof(u64);
		w = (w / n_dicts) * n_dicts;
		if (w == 0)
			errx(1, "RAM budget too small: %" PRIu64 " bytes/node cannot hold one slot per shard",
			     nbytes_memory);
		assert(w % n_dicts == 0);
		w_shard = w / n_dicts;
	}
};

}
#endif
