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

/* what one dict thread got for its shard's pages: written once, before it enters its first round */
struct ShardPages {
	u64 nbytes = 0;                        /* the mapping, the slots rounded up to a whole huge page */
	u64 huge = 0;                          /* of it, what the kernel really backs with 2 MB pages */
	int hugetlb_errno = 0;                 /* why MAP_HUGETLB was refused; 0 if it was not */
};

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
	std::vector<ShardPages> pages;         /* per dict thread, for thread 0 to report */
	Atomic<u32> pages_ready{0};            /* dict threads that have published theirs */

	Shared(int n_threads, int n_dicts) : tally(n_threads), pages(n_dicts) {}

	/* a dict thread's shard, once it is mapped and touched: thread 0 reports them all together */
	void publish_pages(int shard, u64 nbytes, u64 huge, int hugetlb_errno)
	{
		pages[shard].nbytes = nbytes;
		pages[shard].huge = huge;
		pages[shard].hugetlb_errno = hugetlb_errno;
		pages_ready.fetch_add(1);          /* the entry is written before the count that publishes it */
	}

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
