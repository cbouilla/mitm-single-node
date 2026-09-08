#ifndef MITM_ROUTER
#define MITM_ROUTER


/*
 * Route points from sender threads to globally numbered receiver threads across MPI ranks.  Procedural
 * API around one RAII handle, the thread's own Router_thread.  This file is the umbrella; the parts below
 * include each other by relative path, and one of them is include/router/placement.hpp -- NOT the engine's
 * include/placement.hpp, which a bare "placement.hpp" from this directory would not reach.
 *
 * Setup, collective over every thread of the team and every node, the same arguments on every thread but
 * `group` (a thread's own):
 *   Router_thread rt = Router_Init(role, group, comm, tag, lossy, opts)
 *       role    ROUTER_SENDER or ROUTER_RECEIVER; thread 0 is the service whatever role says
 *       group   ROUTER_GROUP_AUTO to let the Router form the groups, or this worker's color in 0..G-1
 *
 * The team's shape and this thread's place in it (any thread, once connected):
 *   Router_local_rank / Router_rank               its rank among its role, on this node / on all nodes
 *   Router_local_num_send / _num_recv             senders, receivers on this node
 *   Router_num_send / Router_num_recv             senders, receivers over all nodes
 *   Router_group / Router_num_groups              its group (a sender's or receiver's; -1 service) / this node's G
 *   Router_group_num_receivers / _receiver(i)     the receivers of its group: a producer pairs with one of them
 *   Router_domain / Router_num_domains            its cache domain (-1 when not pinned) / this node's count
 *   Router_cpu / Router_numa_node                 where the kernel put it
 *
 * The round (sender / receiver / service):
 *   Router_Push(a, b, dest, rt)                   sender: the point (a, b) to receiver `dest`
 *   Router_Close(rt)                              sender: nothing more until Reset
 *   Router_Grab(&pts, rt) / Router_Release(rt)    receiver: the next block's points in place, how many; done with it
 *   Router_Pop(&a, &b, rt)                        receiver: one point at a time, on top of the two above
 *   Router_Test_drained(rt)                       receiver: nothing more will come and nothing is held
 *   Router_Progress / Router_Test_quiescent       service: one bounded, non-blocking turn in a loop; then done?
 *   Router_Stats / Router_Reset                   service: the node's tallies; then a fresh round (an MPI_Barrier)
 *
 * The round, per node: senders push then close; receivers grab, read in place and release until drained; the
 * service calls Progress until quiescent; a barrier; Stats, Reset; a barrier.  Lossy: no call ever waits, a
 * point that cannot be routed at once is dropped and counted.  Lossless: a sender may block, no point is ever lost.
 *
 * Placement (Router_Opts): by default the Router pins the threads and forms the groups.  A group is senders
 * and receivers in one cache domain, at most `group_size` cores each, and the receivers and senders are
 * round-robined over the domains into the least-filled group of each; the staging will shard on it, and a
 * consumer pairs a producer with the dict thread of its group.  `pin = false` leaves the
 * CPUs to the caller and groups by worker index; a per-thread `group` (a color) takes the caller's groups, as
 * MPI_Comm_split does.  `cache_level` overrides the level a domain sits in (0: the lowest shared by several
 * cores).  When pinning, the Router refuses a mask too small for the team or shared with a co-hosted rank.
 */


#include "common.hpp"
#include "connect.hpp"
#include "workers.hpp"
#include "service.hpp"

#endif
