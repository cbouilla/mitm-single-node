#ifndef MITM_DIRECT_COMMON
#define MITM_DIRECT_COMMON

#include <mpi.h>
#include <sched.h>
#include <err.h>
#include <atomic>
#include <mutex>
#include <vector>

#include "tools.hpp"
#include "parameters.hpp"
#include "router/router.hpp"

/*
 * The direct scheme's data: its parameters, its counters and the round record they travel in, and the state
 * a rank's threads share.  The two worker rounds are direct_producer.hpp and direct_dict.hpp; the service
 * round, the epilogue, the printing, the problem wrappers and the entry points are direct.hpp.
 */

namespace mitm::direct {

/* the two phases of a direct round: f fills the dictionary, then g probes it */
enum phase {FILL, PROBE};

/* the tag of every Router message on the communicator; nothing else on it carries this tag */
static constexpr int ROUTER_TAG = 1;


/********************************* parameters ********************************/

/* the CPUs of the calling thread's affinity mask: what the team fills when --producers-per-node is 0 */
static int affinity_cpus()
{
	cpu_set_t mask;
	CPU_ZERO(&mask);
	if (sched_getaffinity(0, sizeof(mask), &mask) != 0)
		err(1, "direct: sched_getaffinity");
	return CPU_COUNT(&mask);
}

/*
 * Everything a run needs beyond the Options, derived once before the team and const from then on: the MPI
 * topology, the team's shape, the dictionary's slots, the domain and how many rounds cover it.  Data only.
 */
struct Params : Options {
	int rank;                              /* this node: one rank per node */
	int n_nodes;                           /* MPI ranks */
	int S;                                 /* producers per node: the Router's senders */
	int R;                                 /* dict threads per node: the Router's receivers */
	int n_threads;                         /* 1 + R + S: the OpenMP team */
	int n_producers;                       /* producers over all nodes */
	int n_dicts;                           /* dict threads over all nodes */
	u64 nbytes_memory;                     /* RAM budget per node, for the shards */
	u64 w;                                 /* slots in the whole, distributed dictionary */
	u64 w_shard;                           /* slots per shard */
	int n;                                 /* domain bits: a preimage is n-bit */
	int m;                                 /* range bits: an image is m-bit */
	u64 domain;                            /* 2^n */
	u64 per_round;                         /* entries a round inserts: fill * w, at most the domain */
	u64 n_rounds;                          /* ceil(domain / per_round), capped by --nrounds */
	bool capped;                           /* --nrounds cut the round count: the search is not exhaustive */
	int check_bits;                        /* bits of the mixed key a slot keeps beside the preimage: 63 - n */

	Params(const Options &o, u64 nbytes_memory, int n, int m) : Options(o), nbytes_memory(nbytes_memory), n(n), m(m)
	{
		MPI_Comm_rank(mpi_comm, &rank);
		MPI_Comm_size(mpi_comm, &n_nodes);
		router.verbose = router.verbose && verbose;   /* the Router prints on rank 0 by itself */
		verbose = verbose && (rank == 0);
		if (n > 63)
			errx(1, "direct: a %d-bit preimage leaves an 8-byte slot no room for its occupancy bit (n <= 63)", n);
		if (m > 64)
			errx(1, "direct: %d-bit images do not fit a 64-bit key", m);
		if (dicts_per_node < 1)
			errx(1, "--dicts-per-node %d: at least one dictionary shard per node", dicts_per_node);
		R = dicts_per_node;
		S = producers_per_node;
		if (S == 0)
			S = affinity_cpus() - 1 - R;
		if (S < 1)
			errx(1, "direct: no CPU left for a producer (%d in the mask, 1 service thread, %d dict threads)",
			     affinity_cpus(), R);
		producers_per_node = S;
		n_threads = 1 + R + S;
		n_producers = n_nodes * S;
		n_dicts = n_nodes * R;

		if (nbytes_memory == 0)
			errx(1, "the RAM budget per node must be given (and nonzero)");
		w_shard = (nbytes_memory * n_nodes / sizeof(u64)) / n_dicts;
		w = w_shard * n_dicts;
		if (w_shard == 0)
			errx(1, "RAM budget too small: %" PRIu64 " bytes/node cannot hold one slot per shard", nbytes_memory);
		if (not (fill > 0 && fill <= 0.9))
			errx(1, "--fill %.2f: the dictionary fill ratio must be in (0, 0.9]", fill);
		domain = 1ull << n;
		per_round = fill * (double) w;
		if (per_round == 0)
			errx(1, "RAM budget too small: %.2f * %" PRIu64 " slots holds no entry", fill, w);
		if (per_round > domain)
			per_round = domain;
		n_rounds = (domain + per_round - 1) / per_round;
		capped = (max_versions < n_rounds);
		if (capped)
			n_rounds = max_versions;
		check_bits = 63 - n;
	}
};


/********************************** counters **********************************/

/* every tally of the scheme; one u64[N_COUNTERS] per thread, written by its owner, summed by thread 0 */
enum counter {
	N_EVAL = 0,             /* producers: evaluations of f (FILL) or g (PROBE), one point pushed each */
	N_INSERT,               /* dict threads: entries inserted (FILL) */
	N_PROBE,                /* probes retired (PROBE) */
	N_STEPS,                /* slots visited by inserts and probes: the cost of linear probing */
	N_MATCH,                /* slots whose check bits matched a probe */
	BAD_MATCH,              /* ... but whose preimage does not map to the key: a check-bit false positive */
	N_COLLISIONS,           /* ... and whose preimage does: f(x) == g(y), golden or not */
	N_COUNTERS
};

/* the round record a node contributes to the epilogue's MPI_Allgather: counters, Router stats, its golden pair */
enum record {
	REC_ROUTER = N_COUNTERS,                        /* ROUTER_STATS_SIZE words: the node's Router_Stats */
	REC_FOUND = N_COUNTERS + ROUTER_STATS_SIZE,     /* 1 if the node holds a golden pair */
	REC_X,                                          /* its x */
	REC_Y,                                          /* its y */
	REC_WORDS
};


/********************************** shared state **********************************/

/* one thread's tallies, on a cache line of its own: written by the owner, read by thread 0 */
struct alignas(64) Tally {
	u64 ctr[N_COUNTERS];                   /* indexed by enum counter */
};

/* what a rank's threads share: the tallies, the golden pair a dict thread found, the epilogue's verdict */
struct Shared {
	std::vector<Tally> tally;              /* per OpenMP thread id */
	std::mutex golden_mtx;                 /* serialises set_golden */
	std::atomic<u32> found{0};             /* 1 once golden[] holds this node's pair */
	u64 golden[3] = {};                    /* the round, x and y of the first golden pair found on this node */
	u64 stop = 0;                          /* no next phase; thread 0 writes it, everyone reads it after the barrier */
	u64 phases = 0;                        /* phases run to their end, thread 0's count */
	bool solved = false;                   /* the search found a pair, on this node or another */
	u64 solution[3] = {};                  /* its round, x and y, the same on every node */

	Shared(int n_threads) : tally(n_threads) {}

	/* a dict thread's golden pair; the first one wins */
	void set_golden(u64 round, u64 x, u64 y)
	{
		std::lock_guard<std::mutex> lock(golden_mtx);
		if (found.load(std::memory_order_relaxed))
			return;
		golden[0] = round;
		golden[1] = x;
		golden[2] = y;
		found.store(1, std::memory_order_release);
	}
};

}
#endif
