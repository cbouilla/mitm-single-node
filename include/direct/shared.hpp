#ifndef MITM_DIRECT_SHARED
#define MITM_DIRECT_SHARED

#include <mutex>
#include <vector>

#include "tools.hpp"
#include "direct/params.hpp"

/*
 * What a rank's threads share: each thread's tallies, the golden pair a dict thread may find, and the
 * epilogue's verdict once every node has reported in.
 */

namespace mitm::direct {

/********************************** shared state **********************************/

/* one thread's tallies, on a cache line of its own: written by the owner, read by thread 0 */
struct alignas(64) Tally {
	u64 ctr[N_COUNTERS];                   /* indexed by enum counter */
};

/* what a rank's threads share: the tallies, the golden pair a dict thread found, the epilogue's verdict */
struct Shared {
	std::vector<Tally> tally;              /* per OpenMP thread id */
	std::mutex golden_mtx;                 /* serialises set_golden */
	Atomic<u32> found{0};                   /* 1 once golden[] holds this node's pair */
	u64 golden[2] = {};                    /* x and y of the first golden pair found on this node */
	u64 stop = 0;                          /* no next phase; thread 0 writes it, everyone reads it after the barrier */
	u64 phases = 0;                        /* phases run to their end, thread 0's count */
	bool solved = false;                   /* the search found a pair, on this node or another */
	u64 solution[2] = {};                  /* x and y, the same on every node */

	Shared(int n_threads) : tally(n_threads) {}

	/* a dict thread's golden pair; the first one wins */
	void set_golden(u64 x, u64 y)
	{
		std::lock_guard<std::mutex> lock(golden_mtx);
		if (found)
			return;
		golden[0] = x;
		golden[1] = y;
		found.store_release(1);
	}
};

}
#endif
