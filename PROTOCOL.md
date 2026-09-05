# Communication protocol

One engine core, several *schemes*.  The core -- `include/engine.hpp` (the comm thread and the
round loop), `include/comm.hpp` (the contexts, the point buffers), `include/spsc.hpp`,
`include/controller.hpp`, `include/parameters.hpp` -- moves points from producer threads to the
dictionary shards and runs the rounds; it is templated on a **scheme**, which says what a point
is, what a dict thread does with it, what a producer does, when a round ends and how it is
reported (§7).  Today's scheme is PCS: `include/pcs_common.hpp` (its data and its `Scheme`),
`include/walker.hpp`, `include/inserter.hpp` and `include/pcs.hpp` (its functions and its entry
points).  A heading or a paragraph marked **(PCS)** describes that scheme's use of the core;
everything else is the core and holds for every scheme.  Names in code font are the ones used
there.

## 1. Actors

**Nodes.**  One MPI rank per node on the provided communicator
  (`Options::mpi_comm`). Rank 0 additionally hosts the **controller**, which
  decides when a round ends and prints the statistics.

**Threads.**  Each rank runs one OpenMP team of `n_threads = 1 + dicts_per_node +
producers_per_node` threads, laid out by thread id:

| tid                  | role       | does                                                                    |
|----------------------|------------|-------------------------------------------------------------------------|
| `0`                  | `COMM`     | all MPI traffic; routes points; on rank 0, runs the controller          |
| `1 .. R`             | `DICT`     | owns dictionary shard `tid-1` (the scheme's `Dict`); consumes every point delivered to it |
| `R+1 .. n_threads-1` | `PRODUCER` | the scheme's producer.  PCS: walks trails, ships distinguished points, resolves collision candidates |

(`R == dicts_per_node`.)  MPI is initialised with `MPI_THREAD_FUNNELED`
and **only thread 0 ever calls MPI** (its code is the `CommThread` class).

**One rank per NUMA node.**  The engine assumes it, and rank 0 prints a loud warning
when its own affinity mask spans more than one (`Placement::report()`, in the banner).  Nothing
enforces it and nothing breaks without it; the placement below simply stops meaning
what it says, because a thread group would then straddle two memory domains.

**Thread groups.**  The mask's **cores** are partitioned into **one group per
dict thread**, each group inside a single cache domain and all of them the same size to
within one core.  The domain is the lowest cache level shared by several *cores* -- L3
on both a Xeon (one per socket) and an EPYC (one per CCX), never hardcoded, detected by
`topology_of_cpus()` and overridable with `--cache-level`.  Groups are spread over the
domains in proportion to their cores, and a domain's leftover cores go to its earliest
groups, because group 0 hosts the comm thread as well as its dict thread.  **Cores, not
CPUs**: a core's SMT siblings have to land in the same group, or the sibling of one
group's dict thread would belong to another group's producer.  Two configurations do not fit
and rank 0 warns: more domains than dict threads, where a group covers a run of whole
domains instead of sitting inside one; and a group left with no CPU for a producer, which
then borrows one from the fullest group.  Fewer producers than dict threads is refused
outright when the scheme asks for a producer in every group (`producer_per_dict`); PCS does,
because every collision queue needs a producer to drain it (§3.2, §4.5).

**Placement, and what SMT is for.**  All of this lives in `placement.hpp`, in the
`Placement` that `Parameters` holds; the engine reads `thread_cpu` and `thread_group`
from it and nothing else.  Every thread takes the **emptiest core of its group**, the service threads first: the comm thread, then dict thread `j` in group `j`,
then the producers spread over the groups, emptiest group first.  Two properties come out
of that one rule.  Each service thread lands on a core of its own while cores last --
two threads that spend their time waiting, one on queues and one on random DRAM, would
only slow each other down on one core.  And the SMT siblings they leave are the first
CPUs the producers fill, which is what we want: a producer is a dependent chain of
vectorised evaluations, so it fills exactly the issue slots a waiting thread leaves
idle.  Measured (`PROBLEM.md` §3): idling those siblings is worth ~10% to the service
threads, while SMT is worth +66% to the producers -- so they are paired, never idled, and
never given to a second service thread.  When a group has fewer cores than service
threads the rule degrades to sharing and rank 0 warns.  Threads beyond the mask's CPUs
stay unpinned and are spread over the groups by index, and the rank warns.  `--no-bind`
(`Options::bind_threads == false`) leaves every thread unpinned and assigns producer `s`
to group `s mod dicts_per_node`.  Once pinned, each thread records the CPU and NUMA
node the kernel reports (`ThreadContext::cpu`, `numa_node`) and its group
(`ThreadContext::group`).

**A producer's group IS the dict thread it works for.**  That is the whole point of the
partition: one dict thread per group means the producer-to-dict-thread map needs no table.
Worker threads communicate with thread 0 through per-thread
Single-Producer-Single-Consumer (SPSC) queues, one `state` atomic each and their
plain-`u64` tallies (all three in their `ThreadContext`; thread 0's own `state` is the
phase of its round loop, §3.3), and with each other only through whatever the scheme hangs on
the `SharedContext` -- for PCS, the collision queue of their own group
(`SharedContext::scheme.coll_q[group]`, §3.2).  Nothing else is shared between workers: PCS's
HyperLogLog registers are one plain array per producer (§3.4).

```
   one node (one MPI rank), one group per dict thread
   +------------------------------------------------------------------+
   |  group 0:  producer --SPSC--+                                    |
   |            producer --SPSC--+-->  comm  --SPSC-->  dict thread 0 |
   |               ^                    |                   |         |
   |               +---- coll_q[0] <----|-------------------+         |
   |                                    |                             |
   |  group 1:  producer --SPSC--+------+    --SPSC-->  dict thread 1 |
   |               ^                    |                   |         |
   |               +---- coll_q[1] <----|-------------------+         |
   +------------------------------------|-----------------------------+
                                        |
          TAG_POINTS      bulk DPs + sentinels  <-->  every rank (self included)
          TAG_END_ROUND   rank 0  -->  every node (self included)
          TAG_REPORT      every node  -->  rank 0
          TAG_SOLUTION    every node  -->  rank 0
          collectives     round header, statistics, answer
```

## 2. Messages between nodes

All messages are arrays of `MPI_UINT64_T` on `mpi_comm`.  There are four tags, one per
message kind: `TAG_POINTS` for the data, and `TAG_END_ROUND`, `TAG_REPORT`,
`TAG_SOLUTION` for the control channel, plus a few collectives.  A message is told
apart by its tag, never by its length.

### 2.1 `TAG_POINTS`: points, node to node

**Bulk point message.**  `2k` words, `1 <= k <= buffer_capacity` (default 1500 points, so at
most 24 000 bytes).  Each pair is one `Point` (`parameters.hpp`):

| word | field | meaning |
|---|---|---|
| `2t+0` | `key` | the routing key, in **full**: it is routed on and never recomputed |
| `2t+1` | `val` | one word of payload, the scheme's |

**(PCS) A point is a distinguished point**: `key` is the endpoint (the distinguished point
itself), `val` the chain index `j` in the low `jbits` and the trail length in the `lenbits`
above it.

**(PCS) Why two words and not three.**  A 64-byte write-combining line holds four points
instead of two, which is worth 2x on the producer-side staging a large node needs
(`PROBLEM.md` §9.1) and cuts the message, both SPSC rings and the double buffers by a
third.  The length is what has to give: it gets `lenbits = 64 - jbits` bits by default
(`--dp-len-bits` narrows it, which is how the paths below are exercised at test scale),
it is shipped as `min(len, len_sat)` where `len_sat = 2^lenbits - 1`, and `len_sat` on
the wire therefore means **"at least that long, exact length unknown"**.  A trail always
has `len >= 1`, so a packed word is never zero.  Recovering an unknown length costs the
producer a re-walk of that trail, §3.2.

Sent with `MPI_Isend` from `OutBuffers`, which keeps one double buffer per
destination: `ready[dst]` accumulates, `outgoing[dst]` is in flight.  When `ready`
is full and the previous send to that node has not completed, the point is
**dropped** (the comm thread tallies it, `DROP_OUT`), never blocked on.  A buffer is
never sent empty.  An `Isend` is tested only when its slot is wanted (`rotate()`) and
at the end of the round (`flush_poll()`); MPI progress is per process, so the receive
tests of every loop turn drive the sends as well, and nothing polls them for their own
sake.

Received into a pool of `n_in_buffers` (default 8) always-posted `MPI_Irecv`s with
`MPI_ANY_SOURCE` (`CommThread::in_req`, one buffer `in_data[k]` each).  The comm
thread's `poll_incoming()` runs `MPI_Testsome` over the pool, scatters each completed
buffer to the local dict queues, and reposts the receive.

**End-of-round sentinel.**  A **zero-length** `TAG_POINTS` message, one to every node
(itself included), sent with `MPI_Bsend` by `OutBuffers::send_sentinels()`.  Meaning:
"I have sent you every point of this round."  A receiver counts them in
`CommThread::n_sentinels`; the round's incoming traffic is over when the count reaches
`n_nodes`.  It is sent only after every data `Isend` has completed locally (step 3 of
the drain, §4.5), and MPI messages between one pair of ranks are non-overtaking, so
the sentinel is matched after the data it closes.

**Routing.**  The dictionary is sharded over all `n_dicts = n_nodes *
dicts_per_node` dict threads; a point lives in shard `key % n_dicts`.  The
three divisions are split over the three stages that need them:

| stage | computes | where |
|---|---|---|
| sending comm thread | `node = key % n_nodes` | `route_producer_queues()` |
| receiving comm thread | `slot = (key / n_nodes) % dicts_per_node` | `deliver_local()` |
| dict thread | `key / n_dicts`, the part of the key the shard indexes on | the scheme's dict thread (PCS: `inserter_thread()`) |

A point whose `node` is the local rank never touches MPI: `deliver_local()` pushes it
straight into the dict queue.

### 2.2 The control channel: `TAG_END_ROUND`, `TAG_REPORT`, `TAG_SOLUTION`

Three tags, one per message kind, and one always-posted `MPI_Irecv` per tag on the
rank that consumes it:

| tag | from | to | words | posted by |
|---|---|---|---|---|
| `TAG_END_ROUND` | rank 0 | every node | 0 | `CommThread::req_end_round`, on every rank, with source `0` (no wildcard) |
| `TAG_REPORT` | every node | rank 0 | `N_COUNTERS` | `Controller`, rank 0 only, `MPI_ANY_SOURCE` |
| `TAG_SOLUTION` | every node | rank 0 | `SOL_NWORDS` | `Controller`, rank 0 only, `MPI_ANY_SOURCE` |

The report and solution receives are serviced by the
controller itself (`Controller::service()`); the node's comm code only ever sees the
end-of-round signal.  Rank 0 talks to itself through MPI like any other node.  Every
send is `MPI_Bsend`: it completes locally, so the comm thread never blocks and there is
no request to track.  The attached buffer is per process -- `run()` attaches it before
anything else of the engine exists and detaches it after the last round (§4.1, §4.7) --
and holds `2 * n_nodes + bsend_slack` report-sized slots.  The bounded part is an end-of-round
signal to every node (rank 0) and one sentinel per node, both once per round; our own
reports and solution have no hard bound (reports are one-way and nothing throttles them
but their pacing rule), and `bsend_slack` (default 8) covers them: a report's slot is
reclaimed by the sender's own progress, and at least 256 turns of the comm loop, each
with several `MPI_Test`s, separate two reports.  A full buffer is a fatal
`MPI_ERR_BUFFER`, never a hang.

The three tags are independent channels: MPI orders only the messages that match the
same receive, so a solution may be matched before or after a report the same node sent
earlier.  Nothing depends on that order: both only feed the controller's decision to
close the round, which is taken once per `service()` call after both receives have
been drained.

**Progress report** (node -> rank 0, `TAG_REPORT`, `Scheme::N_COUNTERS` words).  The
node's counters -- every thread's `ctr[N_COUNTERS]` summed (§3.4), in the order of the
scheme's `enum counter` -- as **deltas since the node's previous report**; the
controller sums them into `reported[]`.  The core reads one of them, `Scheme::PACING`, to
pace the reports by volume; the scheme reads them all, in `round_complete()` (which closes
the round) and `display()`.  One layout serves the per-thread tallies, the report and the
end-of-round reduction alike.

**(PCS) The counters** (`pcs_common.hpp`, 17 of them; `PACING` is `N_DP`):

| index | meaning | tallied by |
|---|---|---|
| `N_EVAL` | evaluations of the mixing function, walking or resolving | producers |
| `N_DP` | distinguished points found | producers |
| `N_POINTS_TRAILS` | sum of the lengths of the trails that reached a DP | producers |
| `N_COLLISIONS` | collisions located | producers |
| `COLLIDING_LEN_MIN` | sum of the shorter trail length of each colliding pair | producers |
| `COLLIDING_LEN_MAX` | ... and of the longer one | producers |
| `N_MEASURE` | trails re-walked because their length had saturated (§3.2) | producers |
| `BAD_DP` | trail gave up before a distinguished point, walking or re-walked to be measured (§3.2) | producers |
| `BAD_COLLISION` | the two trails "collide" on the same value | producers |
| `BAD_WALK_ROBINHOOD` | one trail is a suffix of the other | producers |
| `BAD_WALK_NONCOLLIDING` | dictionary false positive: the trails never meet | producers |
| `DROP_PRODUCERQ` | points dropped: producer queue to comm full | producers |
| `N_PROBE` | dictionary probes retired | dict threads |
| `BAD_PROBE` | probe missed: slot empty, or held a different key | dict threads |
| `DROP_COLL` | candidates dropped: collision queue full | dict threads |
| `DROP_OUT` | points dropped: outgoing MPI buffer busy | comm thread |
| `DROP_DICTQ` | points dropped: dict queue full | comm thread |

Pacing: the comm thread considers reporting every 256 turns of its poll loop, and
sends a report if either `ping_delay` seconds (default 0.1) have elapsed since the
last one, or the `PACING` counter has grown by `report_points` since the last one.  The
scheme's `Params` sets `report_points` (PCS: `points_per_version / (reports_per_round *
n_nodes)` DPs, so that a short round cannot overshoot `beta * w`).  A report
is **one-way**: nothing is answered, and a node keeps reporting through its drain,
until its comm thread goes quiescent (§4.5): no DP is produced past step 1, so those
reports come one per `ping_delay` and carry the collisions resolved and the last
points probed while the round winds down.

**Solution** (node -> rank 0, `TAG_SOLUTION`, `SOL_NWORDS == 3` words).  `SOL_I`,
`SOL_X0`, `SOL_X1`: the mixing-function index and the two colliding points, exactly
the layout of `SharedContext::golden`.  Sent once per node (`golden_sent`) as soon as the
comm thread sees `SharedContext::golden_found`, in whatever phase of the round.  It is
**one-way**.  The controller keeps the first one and raises `stop`.

**End of round** (rank 0 -> every node, `TAG_END_ROUND`, **zero-length**).  Meaning:
"stop producing points and drain: the round is over".  Sent by
`Controller::service()` to every node, rank 0 included, as soon as it sees `stop`
raised or the scheme's `round_complete(reported)` true (PCS: the summed `N_DP` of the
round's reports at or above `points_per_version`).
It is the only message the controller ever sends, it is unsolicited, and it is sent
**exactly once per round** (`round_closed`): the receive that matches it names source
0 and this tag, so successive signals are non-overtaking and a node consumes exactly
one per round -- it is its only way out of the steady state -- and a second one would
be consumed in the *next* round and end it at once.  Closing the round is also where
the search advances: the controller raises `stop` if the round just closed was the
`max_versions`-th (`nround` counts the rounds started), so that `stop` is final by the time
the next round header is broadcast (§4.2).  The controller only knows what has been
reported, so a node that has not reported cannot push `reported[]` over the threshold; that
is why reports are also paced by volume.

### 2.3 Collectives

All on `mpi_comm`, all issued by thread 0, all `MPI_UINT64_T` unless stated.

| when | call | payload | root |
|---|---|---|---|
| driver `mitm::init`, if the seed was not given | `MPI_Bcast` | 1 word: PRNG seed (this is in `examples/driver.hpp`, not the library) | 0 |
| `run()`, before the team exists | `MPI_Bcast` | 3 words: `a`, `b`, `wrapper.self_test(a, b)` (PCS: `mixf`); every rank asserts it computes the same value | 0 |
| start of every round | `MPI_Bcast` | the scheme's `Header` as `u64` words, then `stop` (PCS: `i`, the function version, `root_seed`, `stop`: 3 words) | 0 |
| end of every round | `MPI_Reduce`, `MPI_SUM` | `Scheme::N_COUNTERS` words: the node's `ctr` arrays summed over its threads, in `enum counter` order (the layout of a report, totals instead of deltas) | 0 |
| end of every round | the scheme's `RoundStats::reduce()` | (PCS) `MPI_Reduce`, `MPI_MAX` of `HLL_REGISTERS == 65 536` `MPI_UINT8_T`: the node's HyperLogLog of this round's collisions (its producers' arrays merged once every one of them is quiescent) | 0 |
| end of every round | `MPI_Gather` | 4 words per node: `found`, `i`, `x0`, `x1`, the node's golden slot (§4.6) | 0 |
| after the last round | `MPI_Bcast` | 4 words: `found`, `i`, `x0`, `x1` | 0 |

These collectives are the **only blocking MPI calls** the engine makes after
startup.  The per-round ones are reached only after every node has finished its drain
(§4.5), so no point-to-point traffic of the round can still be pending when they run.

## 3. Channels inside a node

| channel | producer | consumer | type | capacity (`Options`) | when full |
|---|---|---|---|---|---|
| `shared.ctx[R+1+w]->q` | producer `w` | comm | `SPSCQueue` of `DP` | `producer_queue_capacity` (1024) | producer drops the DP, tallies `DROP_PRODUCERQ` |
| `shared.ctx[1+r]->q` | comm | dict thread `r` | `SPSCQueue` of `DP`, filled from a per-shard bucket (§3.1) | `dict_queue_capacity` (4096) | comm drops the points of the run that did not fit, tallies `DROP_DICTQ` |
| (PCS) `shared.scheme.coll_q[r]` | dict thread `r` | the producers of group `r` | `CollisionQueue` (mutex, bounded), transferred in runs both ways | `coll_queue_capacity` (8192), **per dict thread** | dict thread drops what did not fit, tallies `DROP_COLL` |
| `ThreadContext::state` | comm and the thread itself (thread 0's: itself alone) | the other side | `atomic<int>` | | see §3.3 |
| `ThreadContext::ctr` | the thread (the comm thread's own drops in `ctx[0]`) | comm | plain `u64[N_COUNTERS]`, no atomics, see §3.4 | | |
| `SharedContext::header`, `stop` | comm (thread 0) | everyone | the scheme's `Header` and a plain `u64`, published by an OpenMP barrier | | |
| (PCS) `ThreadContext::scheme.hll` | the producer itself | comm, after the round | plain `u8[HLL_REGISTERS]`, no atomics, see §3.4 | | never full |
| `SharedContext::shards[r]` | dict thread `r` | dict thread `r` | the scheme's `Dict`, built, fed and flushed by its dict thread alone.  PCS: `PcsDict` | `w_shard` slots | a dictionary.  PCS: a slot is overwritten by a trail at least as long (`pop_insert`), a trail of unknown length counting as the longest |
| `SharedContext::golden` | any producer | comm | mutex + `atomic<bool>` flag | | first one wins |

Where state lives, per rank: private to a worker thread and never touched by the comm
thread, a local of that thread's function; specific to one thread but also touched by
the comm thread (`role`, `state`, `ctr`, the SPSC queue), its `ThreadContext`; common
to all threads, the `SharedContext` (the round header, the `ctx`, `shards` and `coll_q`
tables, the golden pair), one per rank, a local of `run()`; private to the comm thread
(MPI buffers, requests, the controller), a member of `CommThread`.  What a scheme adds to a
context is its `scheme` member: `ThreadContext::scheme` is the scheme's `ThreadStats` (PCS:
a producer's HyperLogLog), `SharedContext::scheme` its `Shared` (PCS: the collision queues).
Every counter, the comm thread's own included, is in a `ThreadContext`, so that the
node's tallies are one loop with no special case (§3.4); the shards (and PCS's collision
queues) are in the `SharedContext` although only their dict thread produces into them,
because zeroing a shard may one day be collective and because a queue's consumers are a
whole group.

Each SPSC queue belongs to the `ThreadContext` of its worker side (the producer, or
the dict thread), which builds it once that thread is pinned (§4.1).  SPSC capacities are rounded up to a
power of two.  Every queue is **non-blocking on both sides**: nobody ever waits for
room or for data, they poll and move on.

### 3.1 The SPSC queues

Lamport's single-producer/single-consumer ring: the producer owns `tail`, the
consumer owns `head`, each side caches the other's index and reads the shared one
(with acquire) only when its cache says the queue is full or empty.  Publication is a
release store of the owned index.  The comm thread drains producer queues in batches of
64 (`pop_bulk`); dict threads drain theirs in batches of 64 as well.  `head` and `tail`
are public so that the comm thread can test `head == tail` from outside at the end
of a round (§4.5, step 2).

A dict queue is filled in **runs**, not point by point.  The comm thread stages
each local point in a bucket of `CommThread::BUCKET_CAP` (128) per shard and calls
`push_bulk` -- one release store, contiguous slots -- when a bucket fills or when the
call that filled it returns.  One sweep of the producer queues scatters over every shard
of the node, so a push per point would make each ring's lines ping-pong between the
producer and the consumer once per point; a run of `k` costs one store for `k` points.
`push_bulk` pushes a prefix and reports how many fit, so a full queue still drops the
tail of the run rather than blocking.  **The buckets never outlive the call that filled
them**: `route_producer_queues()` and `poll_incoming()` each flush every bucket before
returning, so no point is held across a turn of the comm loop and the drain (§4.5) is
unchanged -- the phase tests of steps 2 and 6 see exactly what they saw before.

### 3.2 The collision queues (PCS)

`CollisionCandidate` is what a dict thread hands to the producers of its own group, through
its own `SharedContext::scheme.coll_q[r]`, on a dictionary hit:

| field | meaning |
|---|---|
| `i` | function version; `service_collision` asserts it equals the current round's |
| `seed0`, `len0_maybe` | the incoming point's chain index and trail length; `len0_maybe == 0` means the length had saturated **on the wire** (§2.1) |
| `end` | the shared endpoint, in full: what a re-walked trail must reach |
| `seed1`, `len1_maybe` | the point already in the slot; `len1_maybe == 0` means the stored length had saturated (8 bits) |

**Either length can be unknown, or both.**  A slot holds 8 bits of length, the wire
holds `lenbits`, and `0` means "unknown" on both sides -- free as a marker because a
trail's length is at least 1.  A length that saturated on the wire is stored as unknown
whatever the 8-bit field could have held, so the two are consistent: the dictionary never
reports an exact length it does not have.  The producer recovers an unknown length by
re-walking that trail to its distinguished point, which must be `end` and must be reached
within `dp_max_it` steps; `N_MEASURE` counts the re-walks, `BAD_DP` and
`BAD_WALK_NONCOLLIDING` the two ways one fails.

**One queue per dict thread, and candidates cross it in runs.**  Collisions are not rare
events: at `beta = 8` a fair fraction of every probe hits, so the queue carries a fixed
share of the traffic the shards do (`PROBLEM.md` §8.5, §9).  One queue for the whole
node, pushed by every dict thread and polled by every producer, is what that share cannot
survive -- so there is one per dict thread, its consumers are the producers of that
dict thread's group, and both sides move up to 64 candidates per lock acquisition (a
dict thread fills a private run and hands it over as soon as it has nothing left to probe;
a producer takes a run into its resolver and works from that).  The mutex stays: what
made the single queue costly was that there was one of it, and that its `count` line was
written millions of times a second, so every producer's poll missed.  Batched and split,
that line is written thousands of times a second and a poll finds it valid.  The
emptiness probe a producer runs before taking the lock is a relaxed atomic read, so an idle
queue still costs no lock traffic at all.

**How a producer retires them.**  A round spends about `2/beta` of its evaluations
locating collisions, and resolving one point at a time makes each of those cost `vlen`
times what a walked one does -- most of a producer's time, once the dictionary fills.  So
a producer with `vlen > 1` keeps `vlen/2` candidates in flight in a `VecResolver` and
steps them all with a single `vmixf`.  A candidate walks its two chains through three
phases, one lane per chain that is moving:

| phase | what it does |
|---|---|
| `MEASURE` | entered when either length is unknown: re-walk that chain -- or both at once -- to its distinguished point to learn its length.  A chain gives up after `dp_max_it` steps (`BAD_DP`), and abandons the candidate if it lands on a point other than `end` (`BAD_WALK_NONCOLLIDING`).  A chain that arrives first simply stops being committed while the other finishes |
| `ALIGN` | rewind both chains and step the longer one until both are the same distance from their shared endpoint; the other waits at its start |
| `MARCH` | step both and compare, until they meet (the collision) or the shorter trail runs out (`BAD_WALK_NONCOLLIDING`).  Equal points on entry mean one trail is a suffix of the other (`BAD_WALK_ROBINHOOD`) |

**A busy slot owns two lanes from `fill()` to `release()`**, whatever phase it is in, so
`n_free == 2 * (nslots - n_busy)` at all times and an empty slot always finds its pair:
with `nslots == vlen/2` the pool cannot run dry.  A chain that is not moving is still
stepped by `vmixf`, its lane simply not committed -- which costs nothing, since a step
is a full-width `vmixf` however many lanes are live.  Every length is exact by the time
`ALIGN` begins, so the march's step budget `min(len0, len1)` is never zero.  Recovering
an unknown length by re-walking the trail rather than recording it costs a few
evaluations and bounds what a lane needs to hold.  `N_EVAL` counts the evaluations a one-at-a-time
resolver would have made, not the lanes spent, so it stays comparable across the two
paths; the idle lanes are the price of the batch and about a third of them are idle.

A step costs a whole `vmixf` however few candidates are in flight, so in the steady
state a producer stops as soon as it has nothing more in hand -- its private run empty
*and* its queue empty -- while its batch is not full, and leaves the partial batch
parked for the next chunk.  **Only the drain (§3.3, §4.5 step 7) runs the batch down to
the last candidate**, where every one of them must be retired whatever it costs.  A
parked candidate holds nothing but its own two chain indices, so parking it is free.

A problem with `vlen == 1` has no vector implementation to batch: its producers resolve
one candidate at a time, in `resolve_collision`.  They draw from the same private run,
which lives in the `VecResolver` for both paths.  That path makes the opposite trade on
an unknown length: `walk_recorded()` **records** the trail it re-walks, in one buffer of
`dp_max_it + 1` points owned by the producer thread, so its march reads the recorded chain
and costs one evaluation per step instead of two.  It records the chain whose length is
unknown and steps the other; when both are unknown, `measure_trail()` re-walks the
stepped one first, without recording it, because `walk_recorded()` has to know how far to
step it.  The march is bounded by the shorter trail, exactly as `MARCH` is.

### 3.3 Per-thread wind-down state

Each worker thread carries a `state` that the comm thread *writes to ask* and the
thread itself *writes to answer*; a thread always announces its own completion.

```
producer:     RUNNING --[comm]--> HOLD --[self]--> HELD --[comm]--> DRAIN --[self]--> QUIESCENT
dict thread:  RUNNING -----------------------------------[comm]--> DRAIN --[self]--> QUIESCENT
```

| state | producer | dict thread |
|---|---|---|
| `RUNNING` | produces points (PCS: walks chunks, ships DPs, retires collision candidates between chunks) | consumes everything in its queue (PCS: probes it) |
| `HOLD` | "stop producing points": acknowledges by writing `HELD` | (not used) |
| `HELD` | no new points; keeps retiring collision candidates, batch still parked between chunks | (not used) |
| `DRAIN` | "nothing more will ever come": finishes what it holds (PCS: empties **its own group's** queue and its private run **and runs its resolver batch down to the last candidate**, §3.2), writes `QUIESCENT`, returns | "no point will ever be delivered": empties its queue, writes `QUIESCENT`, returns |
| `QUIESCENT` | done for this round | done for this round |

A producer reads its state **once per chunk** (`chunk_size` iterations of `vmixf`,
default 64, plus whatever candidates it retires before the chunk), which is why
`HELD` exists: until a producer has acknowledged, an empty producer queue is not proof
that it is idle.  A dict thread checks its state **before** re-testing emptiness, so a
point pushed between the two checks is still probed before it goes quiet.  It hands its
private run of candidates over whenever it has nothing left to probe, so by the time it
answers `QUIESCENT` that run is already in the queue and the drain has nothing to flush.
All state
loads are acquire, all stores release.  Every thread resets itself to `RUNNING` at the
top of the next round.

The comm thread's own `state` (`ctx[0]->state`) is the phase of its round loop (§4.3,
§4.5).  It alone reads and writes it (relaxed), and no worker ever sees these values:

```
comm:      RUNNING -> COLLECTING -> FLUSHING -> WAITING -> DRAINING_DICTS -> DRAINING_PRODUCERS -> QUIESCENT
```

| state | the node is | left when (§4.5) |
|---|---|---|
| `RUNNING` | in its steady state | the end-of-round signal has arrived (step 1) |
| `COLLECTING` | routing what the producers, told to `HOLD`, left in their queues | every producer is `HELD` and every producer queue is empty (step 2) |
| `FLUSHING` | sending its partial output buffers | nothing is left to send and every `Isend` has completed; the sentinels go out (steps 3, 4) |
| `WAITING` | taking delivery until every node's sentinel | `n_sentinels == n_nodes`; the dict threads are told to `DRAIN` (steps 5, 6) |
| `DRAINING_DICTS` | waiting for the dict threads to empty their queues | every dict thread is `QUIESCENT`; the producers are told to `DRAIN` (steps 6, 7) |
| `DRAINING_PRODUCERS` | waiting for the producers to empty the collision queues | every producer is `QUIESCENT` (step 7) |
| `QUIESCENT` | done for this round: `comm_round()` runs the epilogue (§4.6) and returns | |

### 3.4 Statistics and the golden pair

- Every tally of a thread lives in its `ThreadContext::ctr`: one `u64[N_COUNTERS]`
  indexed by `enum counter` (`comm.hpp`).  The owner writes it with **no
  synchronisation at all** -- plain increments, no atomics -- and the comm thread's
  own drops go into `ctx[0]->ctr` the same way: no counter lives anywhere else.  The
  comm thread reads them twice.  During the round, `CommThread::snapshot()` sums the
  arrays of every thread after a `#pragma omp flush`, to build a progress report: it
  sees whatever has reached memory, so a count may lag a little (one bumped once per
  chunk can sit in a register until the producer's next release store, or its next turn on
  its group's collision-queue mutex), which
  is fine for a report.  Once every worker is `QUIESCENT` the values are exact -- a
  worker's last increment happens-before its `QUIESCENT` store (release), which the
  comm thread loaded with acquire (§4.5, steps 6 and 7) -- and they are summed across
  threads, reduced across nodes, and zeroed by the comm thread (§4.6).
- A scheme's statistics beyond the counters are its `ThreadStats` (one per thread, in
  `ThreadContext::scheme`, plain, written by its owner alone) and its `RoundStats`: the comm
  thread `collect()`s the threads' into one once every worker is `QUIESCENT` (which zeroes
  them), `reduce()`s it to rank 0, and the controller `fold()`s it into its all-time copy
  (§4.6).  **(PCS)** the HyperLogLog over the round's collisions is **one per producer**,
  `ThreadContext::scheme.hll`: `HLL_REGISTERS == 65 536` plain `u8`, written by its owner and
  nobody else.  On every collision a producer calls `hll_record(hll, x0, x1)`: the register
  is the top 16 bits of the pair's hash, its value the position of the lowest set bit,
  raised to it and never lowered.  No atomics and no CAS -- the estimator only ever wants
  the maximum of each register, and a maximum is associative, so merging the producers'
  arrays gives exactly what one shared array under CAS would have held.  The comm thread
  merges them (`hll_merge`) once every producer is `QUIESCENT`, reduces the result with
  `MPI_MAX` and zeroes them; the producers' release stores of `QUIESCENT` order the
  registers for its acquire loads.  The merge is `HLL_REGISTERS * producers` byte
  comparisons once per round, tens of milliseconds against a round of tens of seconds.
- `SharedContext::set_golden(i, x0, x1)` takes a mutex, keeps the first triple, and
  raises `golden_found` with a release store.  The comm thread polls the flag (acquire)
  every turn of its loop, in every phase of the round, and forwards it as a solution
  message.
  Neither `golden_found` nor `golden_sent` is ever reset: a golden pair ends the search.

## 4. Lifecycle and synchronization steps

### 4.1 Startup (before any round)

1. The driver calls `MPI_Init_thread(MPI_THREAD_FUNNELED)`; `run()` refuses a lower
   level.  Every rank must use the same PRNG seed (the driver broadcasts it).
2. `run()`, on the main thread: the scheme's `Params` -- the core's `Parameters` plus its own
   -- is built on every rank from identical inputs (rank and size come from
   `MPI_Comm_rank/size`).  Its `Placement` member
   (`placement.hpp`, which owns everything in §1) is built first: it loads the hwloc
   topology once, cuts the affinity mask into groups and fills `thread_cpu` and
   `thread_group`, and resolves `producers_per_node` when the user left it at 0.  Then
   the engine's `MPI_Bsend` buffer (§2.2) is attached, as a local of `run()`: MPI has one per
   process, so whatever the caller had attached is detached first --
   `MPI_Buffer_detach` waits for that buffer's pending sends -- and remembered, to be
   put back at the end (§4.7).  Nothing may `MPI_Bsend` before this, and nothing does.
   The test-vector `MPI_Bcast` (§2.3) asserts every rank iterates the same function;
   rank 0 prints the banner (`Scheme::banner()`).  Then what the threads share is built, as
   locals of `run()`: the `SharedContext` -- the round header and golden slot, the `ctx` and
   `shards` tables, sized but not filled, and the scheme's `Shared` (PCS: the `coll_q` table,
   likewise) -- and the `CommThread`, thread 0's object: the main thread *is* thread 0 of
   the team to come.  Its constructor posts every receive of the engine: the
   `n_in_buffers` point receives and the end-of-round receive itself, and on rank 0
   the `Controller` it holds posts the report and solution receives.  Nothing big is
   allocated before the team exists: the MPI buffers are a few tens of kB.
3. The OpenMP team starts.  Each thread pins itself to `thread_cpu[tid]` and asks the
   kernel (`getcpu`) which CPU and NUMA node it is on; a thread that could not be
   pinned, or is not on its CPU, says so.  **OpenMP barrier**: if any thread of the
   rank is misplaced, thread 0 -- the main thread, the only one allowed to call MPI --
   `MPI_Abort`s the run; nothing has been allocated yet.  **Then** each thread builds
   its `ThreadContext` into `shared.ctx[tid]` (with the CPU and NUMA node it measured,
   its group, for a worker its SPSC queue, and the scheme's `ThreadStats` -- PCS: a
   producer's zeroed HyperLogLog), and a dict thread also calls `Scheme::build_dict()`, which
   builds its dictionary shard into `shared.shards[tid-1]` (PCS: a `PcsDict` of `w_shard`
   slots, zero-filled, and its collision queue into `shared.scheme.coll_q[tid-1]`).
   Allocating after pinning is deliberate: under Linux's first-touch policy a page
   lands on the NUMA node of the CPU that first writes it, so this puts every object
   on its owner's node.  **OpenMP barrier** publishes the tables; the comm thread
   reads them from its first round on, and on rank 0 prints the measured layout, one
   line per NUMA node and one per cache domain.

### 4.2 Round start

1. Every thread sets its own `state = RUNNING`.
2. **Thread 0 only**, `CommThread::begin_round()`: rank 0 has the scheme draw the next
   header from the previous one (`Scheme::next_header()`; PCS: a fresh `i` and `root_seed`)
   and hands it `controller.stop`, which the scheme may raise but never lower;
   **`MPI_Bcast` of the header's words plus `stop`** from rank 0.  Every rank stores them
   into `SharedContext::header` and `stop`.  On rank 0, `controller.begin_round()` counts
   the round and resets the per-round tallies and `round_closed`, unless `stop` is set.
3. **OpenMP barrier**: publishes `SharedContext` to every thread.
4. If `stop`, every thread leaves the round loop (§4.7).  Otherwise thread 0 enters
   `CommThread::comm_round()`, dict threads `Scheme::dict_thread()` (PCS: `inserter_thread()`),
   producers `Scheme::producer_thread()` (PCS: `walker_thread()`).

### 4.3 Steady state

**Comm thread**, one turn of its loop.  The loop runs from the round start until the
node is quiescent; its body is the same in every phase (the phase is `ctx[0]->state`,
§3.3), and only the exit test at the end of a turn depends on it:

1. `route_producer_queues()`: pop up to 64 DPs from each producer queue; local points go
   to `deliver_local()`, which stages them per shard (§3.1), remote ones to
   `outbuf.push()`.  Every bucket is flushed before the call returns.
2. `poll_incoming()`: `MPI_Testsome` on the receive pool; scatter completed buffers to
   the dict queues through the same per-shard buckets, count sentinels, repost,
   flush the buckets.
3. `service_control()`: `MPI_Test` the end-of-round receive; if it completed, set
   `round_over` (and repost it for the next round).  On rank 0, then
   `controller.service()`: drain the solution receive (the first solution is kept and
   `stop` is raised) and the report receive (each report's deltas are added to the
   round's tallies), then close the round (§2.2) if it is not closed yet and `stop` is
   raised or the scheme's `round_complete()` says so: raise `stop` if this was round
   `max_versions`, `TAG_END_ROUND` to every node.
4. If `golden_found` and not yet sent: `MPI_Bsend` the solution to rank 0.
5. Every 256 turns: `snapshot()` the node's tallies (§3.4) and, if the pacing rule
   (§2.2) says so, `MPI_Bsend` a progress report with the deltas since the previous one.
6. The phase's exit test and the step it triggers (§4.5).  In `RUNNING`: if
   `round_over`, step 1 of the drain -- the producers are told to `HOLD` and the phase
   becomes `COLLECTING`.

Nothing in this loop blocks.

**(PCS) Producer**, one turn: read `state`; retire up to `coll_per_chunk` candidates from its
own group's queue (0 = until it has nothing in hand and the batch is not full, §3.2);
walk one chunk,
pushing each DP found into its SPSC queue (drop if full)
and restarting the chain, tallying into `ctr` as it goes; add the chunk's evaluations
to `ctr`.

**(PCS) Dict thread**, one turn: pop up to 64 DPs, probe each into its shard
(`PcsDict::pop_insert`), stage a `CollisionCandidate` for every hit in a private run;
hand the run to its own queue once there is nothing left to probe, or as soon as it is
64 long (drop what does not fit); then read `state`; if nothing to do, `cpu_relax()`.

**Controller** (inside `service()` on rank 0): a solution is recorded (first one wins)
and raises `stop`; a progress report adds its deltas; then, once per round, the round
is closed -- `stop` raised if this was round `max_versions`, `TAG_END_ROUND` to every
node -- as soon as `stop` is raised or the scheme's `round_complete()` holds (PCS:
`ndp >= points_per_version`).  The live one-line display (`Scheme::display()`) is
refreshed at most every 0.5 s.

### 4.4 Point flow, end to end (PCS)

```
producer finds a DP
   |  SPSC
   v
comm thread (sender):  node = key % n_nodes
   |                                    \
   | local                               \  remote: OutBuffers, MPI_Isend TAG_POINTS
   v                                      v
deliver_local()   <-----------------  comm thread (receiver): poll_incoming()
   |  slot = (key / n_nodes) % dicts_per_node, SPSC
   v
dict thread r:  key / n_dicts, unpack (j, len) from val, probe the shard
   |  hit: CollisionCandidate, staged in a private run of 64
   v
coll_q[r]  -->  a producer of group r: takes a run, then vlen/2 at a time,
                walks both trails, locates the collision, tests the pair
                |  golden
                v
          SharedContext::set_golden  -->  comm thread  -->  solution to rank 0
```

### 4.5 End of round: the drain

The same loop, phase by phase, from the arrival of the end-of-round signal.  Every
turn still runs the full body of §4.3 -- routing, delivery, the control channel, the
golden pair, the paced reports -- and then the exit test of its phase; passing it
runs the step and moves to the next phase.  The order is what makes the round
airtight: no thread declares itself finished while something can still arrive for it.

| step | phase | passes when | then | guarantee once it passes |
|---|---|---|---|---|
| 1 | `RUNNING` | `round_over` | `state(PRODUCER) := HOLD`; phase `COLLECTING` | producers will stop producing at their next chunk boundary; they keep resolving collisions |
| 2 | `COLLECTING` | every producer is `HELD`, this turn's `route_producer_queues()` moved nothing, and every producer queue has `head == tail` | phase `FLUSHING` | every DP this node produced this round has been routed (delivered locally or handed to `OutBuffers`).  The `HELD` load is acquire and a producer's last push happens-before its `HELD` store, so the queue test that follows sees it; a `HELD` producer never pushes again |
| 3, 4 | `FLUSHING` | `outbuf.flush_poll()`: nothing left to send, every `Isend` completed | `outbuf.send_sentinels()`; phase `WAITING` | every DP bound for another node has left this node -- receiving went on throughout, so the peers we waited on could complete their sends too -- and every node (self included) will learn that we are done sending; the sentinel cannot overtake the data it follows |
| 5, 6 | `WAITING` | `n_sentinels == n_nodes` | `state(DICT) := DRAIN`; phase `DRAINING_DICTS` | every node has finished sending to us, and everything they sent has been scattered to the dict queues.  Every node was sent its end-of-round signal before rank 0 could even enter its drain, so nobody is waiting on rank 0 here; it keeps digesting reports so the round's tallies cover what the others produced meanwhile |
| 6, 7 | `DRAINING_DICTS` | every dict thread is `QUIESCENT` | `state(PRODUCER) := DRAIN`; phase `DRAINING_PRODUCERS` | every DP of the round has been probed and every dict thread's private run has been handed to its queue (§3.3); no new candidate can appear in any of them |
| 7 | `DRAINING_PRODUCERS` | every producer is `QUIESCENT` | phase `QUIESCENT`: the epilogue (§4.6), then `comm_round()` returns | every producer has finished what it held.  PCS: every collision queue is empty, and so is every producer's private run and resolver batch, so every candidate of the round has been resolved; each queue has at least one producer of its own group to drain it, which is why fewer producers than dict threads is refused at startup (§1).  Doing this last is what keeps a golden pair found on the very last candidate from being lost |

Running the full body in every phase is harmless.  From `FLUSHING` on the producers are
held and their queues empty, so the routing pass moves nothing (asserted).  Past
`WAITING` no `TAG_POINTS` message of the round remains, and none of the next round can
exist yet: a node sends round `r+1` points only after the round-`r+1` `Bcast` returns,
which needs rank 0's round-`r` reductions to complete, i.e. every node past its own
drain.  The control channel is serviced in every turn: on rank 0 so that the reports
and solutions of nodes still in their steady state are digested (the round's tallies,
and `stop` for the next round header), and on every rank because its `MPI_Test`s
drive the progress our sentinels and rank 0's signals need to actually leave.  No
node's liveness depends on it: every signal of the round has been sent before rank 0
can enter its own drain.

### 4.6 Epilogue

1. Thread 0, at the end of `comm_round()`, once every worker of the node is
   `QUIESCENT` -- their `ctr` arrays and the scheme's per-thread statistics are then exact
   and nobody writes them (§3.4) -- `CommThread::end_round()`: sum the `n_threads` `ctr`
   arrays (`snapshot()`), `collect()` the scheme's `RoundStats` from the threads (PCS:
   merge the producers' HyperLogLogs into one, `hll_merge`, and zero them), **`MPI_Reduce`
   (SUM) of the `N_COUNTERS` words to rank 0**, then the scheme's `RoundStats::reduce()`
   (PCS: **`MPI_Reduce` (MAX) of the `HLL_REGISTERS` bytes**), then **`MPI_Gather` of every
   node's golden slot** -- `(found, i, x0, x1)`.  The slot is final, every worker being
   quiescent, so a pair found during the drain reaches rank 0 here whatever became of its
   `TAG_SOLUTION` message (§6).  On rank 0, `controller.record()` keeps the first pair
   gathered and raises `stop`, and `controller.end_round()` folds the round into its
   all-time `total[]` and `RoundStats` and has the scheme print the round report
   (`Scheme::round_report()`) -- printing only: the round count and the round's closing
   were settled earlier (§2.2).  Then thread 0 zeroes every thread's `ctr` and its own
   `n_sentinels`; nobody touches any of them again before the next round's barrier
   (§4.2, step 3).
2. **OpenMP barrier**: every thread of the node is back from its round function.
   Nothing depends on it; it is an explicit synchronisation point, kept as such.
3. Each dict thread runs `Scheme::after_round()` on its own shard (PCS: zeroes it,
   `flush()`, so that every round starts from an empty dictionary).  The shard is its alone,
   so this needs no synchronisation with anyone, and it is done before the dict thread
   reaches the next round's barrier.

### 4.7 Termination

When the round header carries `stop`, every thread of every rank leaves the loop
after the barrier of §4.2 step 3 (no round is armed).  Thread 0 then runs
`CommThread::finish()`: on rank 0 it packs the solution (`found, i, x0, x1`) and
has the scheme print the final line (`Scheme::done()`); every rank **`MPI_Bcast`s the
4-word answer**, then cancels
its posted receives (the point receives, then `Controller::shutdown`, a no-op off
rank 0, then the end-of-round receive).  Back in `run()`, the engine's `Bsend`
buffer is detached -- this waits for its last buffered sends to be out -- and the
caller's buffer, if there was one (§4.1), is attached again.  `run()` returns the
answer on every rank.  The driver calls `MPI_Finalize`.

## 5. Guarantees and loss semantics

**What may be lost is the scheme's call** (`Scheme::LOSSLESS`).  **(PCS)** distinguished
points and collision candidates are loss-tolerant: losing one costs the work that produced
it and nothing else, so every DP/candidate channel drops on overflow rather than blocking,
and every drop is tallied and shown in the round report (`DROPPED ... producer-queue /
output-buffer / dict-queue / collision-queue`):

| channel | may drop | counted in |
|---|---|---|
| producer -> comm SPSC | yes | `DROP_PRODUCERQ` |
| comm -> remote node (`OutBuffers`) | yes | `DROP_OUT` |
| comm -> dict thread SPSC | yes | `DROP_DICTQ` |
| (PCS) dict thread `r` -> its group's producers (`coll_q[r]`) | yes | `DROP_COLL` |
| control channel (end-of-round signals, reports, solutions) | **never** | buffered send, buffer sized for the bounded traffic plus slack (§2.2); a full buffer is a fatal error, not a drop |
| sentinels | **never** | same buffer |
| collectives | n/a | |

**What is enforced per round.**

- Every point produced in round `r` is either dropped (and counted) or consumed by the
  round-`r` dictionary, on the shard that owns it, before that shard's dict thread goes
  quiescent (drain steps 2 through 6).  No round-`r` point is ever consumed in round
  `r+1`: nothing is delivered after the last sentinel, and the queues are empty when
  the dict threads stop.
- (PCS) Every collision candidate of round `r` is resolved by a round-`r` producer of the
  producing dict thread's group before the producers go quiescent (drain step 7) -- the ones
  in a producer's private run and the ones parked in its resolver batch included; a producer
  asserts `c.i == round.i` as it starts one.
- (PCS) The dictionary is empty at the start of every round (each dict thread flushes its
  shard in `after_round()`, before the next round's barrier).
- The controller sends exactly one message per round, the end-of-round signal to every
  node, and answers nothing; a round can end no other way, and a node leaves the
  steady state exactly once per round.

**Why it does not deadlock.**

- The comm thread never blocks on point-to-point MPI: every send is `Isend` (tested,
  never waited) or `Bsend` (completes locally); every receive is an always-posted
  `Irecv` that is only ever `Test`ed.
- The loop body never changes: while waiting for its own sends to complete (step 3),
  for the peers' sentinels (step 5) or for its own workers (steps 6, 7), a node keeps
  receiving and keeps servicing the control channel, so two nodes waiting on each
  other both make progress.
- Rank 0 closes the round for every node in one go, and observes its own signal only on
  a later turn, so by the time it enters its own drain every node has been sent its
  signal.
- The only blocking MPI calls (the round-start `Bcast` and the epilogue reductions)
  are reached by a rank only after its drain, i.e. after every node has received the
  end-of-round signal and has confirmed with a sentinel; so every rank reaches them.
- Inside a node, no thread ever waits on a queue: a collision queue's mutex is held
  only for the length of one run's copy, never across anything that can block.  The only
  waits are the comm thread spinning on `state` values that the worker threads set for
  themselves at chunk boundaries, and the two OpenMP barriers per round (§4.2 step 3,
  §4.6 step 2), which every thread reaches unconditionally once its round function
  returns.

## 6. Known windows

One place where the protocol is slightly weaker than the rules above suggest, and one
that used to be.  Both are benign for the search itself; they are listed so nobody
rediscovers them.

- **Closed: a golden pair found late in the drain used to be reported one round late, or
  never.**  The comm thread forwards `golden_found` in every phase (§4.3 step 4), but rank 0
  digests solutions in `service()`, which runs in every phase of its own round and not
  between rounds: a `TAG_SOLUTION` landing after rank 0 has left its drain was matched in
  the *next* round, and if that round was the last permitted one, never.  A final
  `service()` in `finish()` could not close this: MPI does not order a node's `Bsend`
  against its later collective contribution.  The epilogue's `MPI_Gather` of every node's
  golden slot (§4.6) does: the slot is final when every worker is quiescent, and the
  gather is a collective every rank contributes to.  `TAG_SOLUTION` is now only what makes
  the exit *early*, by closing the round in which the pair was found.
- **Reports are one-way, so the round's `reported[]` tallies are approximate.**  A
  node keeps reporting until its comm thread goes quiescent, so what is never
  reported is the delta between its last report and that moment: whatever its last
  chunks still produced, and up to a `ping_delay` of drain-time tallies; on top of
  that a snapshot reads the workers' tallies without synchronisation (§3.4).  Only the live
  display and the decision to close the round rest on them; the reductions of §4.6
  count everything exactly, and the round report prints those.  Further, a report and
  the same node's later sentinel match different receives, so the MPI standard does
  not order them: a report could in principle be matched only after rank 0 has run
  the epilogue and reset `reported[]`, and be credited to the next round, which could
  then close early by at most one report per node.  Open MPI matches messages per
  source in sequence order and `service()` drains the report receive to empty on
  every call, so this does not happen in practice.

## 7. Schemes

The core is templated on a `Scheme` struct -- `mitm::pcs::Scheme` today, declared in
`pcs_common.hpp` and defined in `pcs.hpp` -- which names the types and constants below and
defines the functions the core calls.  Nothing in the core knows what a point's payload or a
dictionary slot means, or when a search is over.

| the scheme's | is | PCS |
|---|---|---|
| `Params` | its parameters, derived from the core's `Parameters` (rank, layout, `w` slots) and setting `report_points` (§2.2) | `jbits`, the length field, theta, `points_per_version` |
| `Header` | the round header's words: trivially copyable, whole `u64`s; `stop` travels beside it (§4.2) | `i`, `root_seed` |
| `Dict` | one shard, built by `build_dict()` | `PcsDict`: direct-addressed, a slot overwritten by a longer trail |
| `enum counter`, `N_COUNTERS`, `PACING` | its tallies -- the report and reduction layout -- and the one that paces the reports | 17 counters, `N_DP` |
| `ThreadStats`, `RoundStats` | statistics beyond the counters: a thread's, and the round's with `collect()`, `reduce()`, `fold()` (§3.4, §4.6) | a producer's HyperLogLog; the merged registers |
| `Shared` | state hung on the `SharedContext` | the collision queues (§3.2) |
| `LOSSLESS`, `DROP_OUT`, `DROP_DICTQ` | the overflow policy, and the comm thread's drop tallies (§5) | drops |
| `producer_per_dict` (to `Parameters`) | whether `Placement` must give every dict thread a producer (§1) | yes |
| `next_header()` | rank 0's next header and `stop`, from the previous ones (§4.2) | a fresh `i` and `root_seed` |
| `build_dict()`, `after_round()` | a dict thread's shard, built once pinned (§4.1); what it does after the epilogue barrier (§4.6) | shard and collision queue; `flush()` |
| `producer_thread()`, `dict_thread()` | the two worker rounds (§3.3) | `walker_thread()`, `inserter_thread()` |
| `round_complete()` | closes the round on the reported tallies (§2.2) | `N_DP >= points_per_version` |
| `banner()`, `display()`, `round_report()`, `done()` | all the printing | |
