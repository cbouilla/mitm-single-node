#ifndef MITM_PARAMETERS
#define MITM_PARAMETERS

#include <mpi.h>

#include "tools.hpp"
#include "router/router.hpp"

namespace mitm {

/*
 * Every user knob of every scheme, all with a working default.  A scheme ignores the other's.  Data
 * only; a scheme derives what it needs from an Options and the RAM budget, once, at the top of its run.
 */
struct Options {
	int dicts_per_node = 1;                   /* dictionary shards per node. Users SHOULD set this themselves */
	bool verbose = true;                      /* print progress information (from rank 0 only) */

	/* --- below this line, the defaults should be just fine */

	/* nodes */
	MPI_Comm mpi_comm = MPI_COMM_WORLD;       /* one rank per node */

	/* thread layout inside a node: the team is 1 service thread + dicts_per_node + producers_per_node */
	int producers_per_node = 0;               /* 0 == fill the inherited affinity mask */

	Router_Opts router;                       /* the Router's own: pinning, groups, blocks, credits, its
	                                             verbosity.  The same on every node; the Router checks */

	/* direct */
	double fill = 0.5;                        /* dictionary fill ratio: entries per round == fill * slots */

	/* PCS */
	double theta = -1;                        /* proportion of distinguished points. -1 == auto */
	double alpha = 2.5;                       /* auto-chosen theta == alpha * sqrt(w/n) */
	double beta = 8;                          /* use each function variant for beta*w DPs */
	int dp_lenbits = 0;                       /* bits of trail length shipped with a DP.  0 == all that fit */
	u64 multiplier = 0x2545f4914f6cdd1dull;   /* to generate starting points */
	u64 max_versions = 0xffffffffffffffffull; /* how many rounds (PCS: function versions) before giving up */
	size_t coll_queue_capacity = 8192;        /* candidates buffered by each dict thread for its group */
	size_t coll_per_chunk = 0;                /* candidates a producer retires per chunk; 0 == until its
	                                             batch is not full and the queue is empty */
	size_t chunk_size = 64;                   /* vmixf iterations between queue / phase checks */

	/* pacing */
	double ping_delay = 0.1;                  /* seconds between two refreshes of the live line */
};

}
#endif
