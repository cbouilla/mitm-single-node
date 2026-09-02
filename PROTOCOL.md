# Communication protocol

Everything below is implemented in `include/engine.hpp` (the node and its
round loop), `include/comm.hpp`(queues, buffers, control channel),
`include/spsc.hpp`, `include/controller.hpp`, `include/walker.hpp` and
`include/inserter.hpp`.  Names in code font are the ones used there.

## 1. Actors

**Nodes.**  One MPI rank per node on the provided communicator
  (`Options::mpi_comm`). Rank 0 additionally hosts the **controller**, which
  decides when a round ends and prints the statistics.

**Threads.**  Each rank runs one OpenMP team of `n_threads = 1 + inserters_per_node +
walkers_per_node` threads, laid out by thread id:

| tid                  | role       | does                                                                    |
|----------------------|------------|-------------------------------------------------------------------------|
| `0`                  | `COMM`     | all MPI traffic; routes points; on rank 0, runs the controller          |
| `1 .. R`             | `INSERTER` | owns dictionary shard `tid-1`; probes every point delivered to it       |
| `R+1 .. n_threads-1` | `WALKER`   | walks trails, ships distinguished points, resolves collision candidates |

(`R == inserters_per_node`.)  MPI is initialised with `MPI_THREAD_FUNNELED`
and **only thread 0 ever calls MPI**; the workers never see an `MPI_` symbol.
Worker threads communicate with thread 0 through per-thread
Single-Producer-Single-Consumer (SPSC) queues and a handful of atomics, and
with each other only through the collision queue.

```
   one node (one MPI rank)
   +-----------------------------------------------------------+
   |  walker 1 --SPSC--+                                       |
   |  walker 2 --SPSC--+-->  comm  --SPSC-->  inserter 1       |
   |    ...            |      |    --SPSC-->  inserter 2       |
   |  walker W --SPSC--+      |               ...              |
   |      ^                   |                 |              |
   |      +--------------- coll_q  <------------+              |
   +--------------------------|--------------------------------+
                              |
          TAG_POINTS   bulk DPs + sentinels   <-->  every rank (self included)
          TAG_CONTROL  reports / assignments  <-->  rank 0
          collectives  round header, statistics, answer
```

## 2. Messages between nodes

All messages are arrays of `MPI_UINT64_T` on `mpi_comm`.  There are two tags,
`TAG_POINTS` and `TAG_CONTROL`, plus a few collectives.

### 2.1 `TAG_POINTS`: distinguished points, node to node

**Bulk DP message.**  `3k` words, `1 <= k <= buffer_capacity` (default 1500 DPs, so at
most 36 000 bytes).  Each triple is one `DP`:

| word | field | meaning |
|---|---|---|
| `3t+0` | `seed` | chain index `j` the trail was started from |
| `3t+1` | `x` | the **full** endpoint (the distinguished point itself) |
| `3t+2` | `len` | trail length |

Sent with `MPI_Isend` from `OutBuffers`, which keeps one double buffer per
destination: `ready[dst]` accumulates, `outgoing[dst]` is in flight.  When `ready`
is full and the previous send to that node has not completed, the point is
**dropped** (`OutBuffers::n_dropped`), never blocked on.  A buffer is never sent empty.

Received into a pool of `n_in_buffers` (default 8) always-posted `MPI_Irecv`s with
`MPI_ANY_SOURCE` (`InBuffers`).  The comm thread's `poll_incoming()` runs
`MPI_Testsome` over the pool, scatters each completed buffer to the local inserter
queues, and reposts the receive.

**End-of-round sentinel.**  A **zero-length** `TAG_POINTS` message, one to every node
(itself included), sent with `MPI_Bsend` by `OutBuffers::send_sentinels()`.  Meaning:
"I have sent you every point of this round."  A receiver counts them in
`InBuffers::n_sentinels`; the round's incoming traffic is over when the count reaches
`n_nodes`.  It is sent only after every data `Isend` has completed locally (step 3 of
the drain, §4.5), and MPI messages between one pair of ranks are non-overtaking, so
the sentinel is matched after the data it closes.

**Routing.**  The dictionary is sharded over all `n_inserters = n_nodes *
inserters_per_node` inserter threads; a point lives in shard `x % n_inserters`.  The
three divisions are split over the three stages that need them:

| stage | computes | where |
|---|---|---|
| sending comm thread | `node = x % n_nodes` | `route_walker_queues()` |
| receiving comm thread | `slot = (x / n_nodes) % inserters_per_node` | `deliver_local()` |
| inserter thread | `key = x / n_inserters` | `inserter_thread()` |

A point whose `node` is the local rank never touches MPI: `deliver_local()` pushes it
straight into the inserter queue.

### 2.2 `TAG_CONTROL`: reports and assignments

One always-posted `MPI_Irecv(buf, REP_NWORDS, MPI_ANY_SOURCE, TAG_CONTROL)` per rank
carries the whole control channel in both directions (`ControlChannel`).  Messages
are told apart by **length**: 1 word is an assignment, `REP_NWORDS` words is a report.
Rank 0 talks to itself through MPI like any other node.  Every send is `MPI_Bsend`:
it completes locally, so the comm thread never blocks and there is no request to
track.  The attached buffer holds `3 * n_nodes + bsend_slack` report-sized slots
(worst case: an assignment to every node, our own report, one sentinel per node).

**Progress report** (node -> rank 0, `REP_NWORDS == 11` words).  All counts are
**deltas since the node's previous report**; the controller sums them.

| index | field | source |
|---|---|---|
| `REP_NDP` | distinguished points found | walkers' `n_dp` |
| `REP_NEVAL` | (reserved, always 0) | |
| `REP_DROP_WALKERQ` | points dropped: walker queue to comm full | walkers' `n_drop_walkerq` |
| `REP_DROP_OUT` | points dropped: outgoing MPI buffer busy | `OutBuffers::n_dropped` |
| `REP_DROP_INSERTERQ` | points dropped: inserter queue full | `PcsNode::n_drop_inserterq` |
| `REP_DROP_COLL` | candidates dropped: collision queue full | inserters' `n_drop_coll` |
| `REP_GOLDEN` | 0 | |
| `REP_I`, `REP_X0`, `REP_X1` | 0 | |
| `REP_NPROBE` | dictionary probes retired | inserters' `n_probe` |

Pacing: the comm thread considers reporting every 256 turns of its poll loop, and
sends a report if either `ping_delay` seconds (default 0.1) have elapsed since the
last one, or the node has found `points_per_version / (reports_per_round * n_nodes)`
new DPs since the last one (so a short round cannot overshoot `beta * w`).  **At most
one progress report is outstanding per node**: none is sent while
`awaiting_assignment` is set, and it is cleared by the reply.

**Golden report** (node -> rank 0, `REP_NWORDS` words).  `REP_GOLDEN = 1`, and
`REP_I`, `REP_X0`, `REP_X1` carry the mixing-function index and the two colliding
points; every other field is 0.  Sent once per node (`golden_sent`) as soon as the
comm thread sees `RoundState::golden_found`.  It is **one-way**: no reply, and it does
not touch `awaiting_assignment`.

**Assignment** (rank 0 -> node, 1 word).  The reply to every progress report:

| value | meaning |
|---|---|
| `KEEP_GOING` (0) | carry on with this version of the function |
| `NEW_VERSION` (1) | stop producing points and drain: the round is over |

The controller replies `NEW_VERSION` once the summed `REP_NDP` of the round has
reached `points_per_version`, or once `stop` has been raised by a golden report.  It
decrements `n_active_nodes` per `NEW_VERSION` it hands out, and it never sends
anything unsolicited: a node that has not reported cannot be told to stop, which is
why reports are also paced by DP count.

### 2.3 Collectives

All on `mpi_comm`, all issued by thread 0, all `MPI_UINT64_T` unless stated.

| when | call | payload | root |
|---|---|---|---|
| driver `mitm::init`, if the seed was not given | `MPI_Bcast` | 1 word: PRNG seed (this is in `examples/driver.hpp`, not the library) | 0 |
| `run_engine()`, before the node exists | `MPI_Bcast` | 3 words: `x`, `y`, `mixf(x, y)`; every rank asserts it computes the same value | 0 |
| start of every round | `MPI_Bcast` | 3 words: `i` (function version), `root_seed`, `stop` | 0 |
| end of every round | `MPI_Reduce`, `MPI_SUM` | `ST_NWORDS == 10` words: the round's merged `Counters` (see `round_stat`) | 0 |
| end of every round | `MPI_Reduce`, `MPI_MAX` | 65 536 `MPI_UINT8_T`: HyperLogLog registers of this round's collisions | 0 |
| after the last round | `MPI_Bcast` | 4 words: `found`, `i`, `x0`, `x1` | 0 |

These collectives are the **only blocking MPI calls** the engine makes after
startup.  The per-round ones are reached only after every node has finished its drain
(§4.5), so no point-to-point traffic of the round can still be pending when they run.

## 3. Channels inside a node

| channel | producer | consumer | type | capacity (`Options`) | when full |
|---|---|---|---|---|---|
| `walker_q[w]` | walker `w` | comm | `SPSCQueue` of `DP` | `walker_queue_capacity` (1024) | walker drops the DP, `n_drop_walkerq++` |
| `inserter_q[r]` | comm | inserter `r` | `SPSCQueue` of `DP` | `inserter_queue_capacity` (4096) | comm drops the DP, `n_drop_inserterq++` |
| `coll_q` | every inserter | every walker | `CollisionQueue` (mutex, bounded) | `coll_queue_capacity` (8192) | inserter drops the candidate, `n_drop_coll++` |
| `ThreadContext::state` | comm and the thread itself | the other side | `atomic<int>` | | see §3.3 |
| `ThreadContext::n_dp`, `n_probe`, `n_drop_*` | the thread | comm | relaxed `atomic<u64>` | | |
| `ThreadContext::ctr` | the thread | comm, after the round | plain `Counters` | | |
| `RoundState::i`, `root_seed`, `stop` | comm (OpenMP master) | everyone | plain `u64`, published by an OpenMP barrier | | |
| `RoundState::golden` | any walker | comm | mutex + `atomic<bool>` flag | | first one wins |

SPSC capacities are rounded up to a power of two.  Every queue is **non-blocking on
both sides**: nobody ever waits for room or for data, they poll and move on.

### 3.1 The SPSC queues

Lamport's single-producer/single-consumer ring: the producer owns `tail`, the
consumer owns `head`, each side caches the other's index and reads the shared one
(with acquire) only when its cache says the queue is full or empty.  Publication is a
release store of the owned index.  The comm thread drains walker queues in batches of
64 (`pop_bulk`); inserters drain theirs in batches of 64 as well.  `head` and `tail`
are public so that the comm thread can test `head == tail` from outside at the end
of a round (§4.5, step 2).

### 3.2 The collision queue

`CollisionCandidate` is what an inserter hands to the walkers on a dictionary hit:

| field | meaning |
|---|---|
| `i` | function version; `service_collision` asserts it equals the current round's |
| `seed0`, `len0` | the incoming point's chain index and trail length |
| `end` | the shared endpoint, as the dictionary key `x / n_inserters` |
| `seed1`, `len1_maybe` | the point already in the slot; `len1_maybe == 0` means the stored length had saturated (8 bits), so the walker re-walks that trail (`walk_nolen1`) |

Pushes and pops take the mutex; the emptiness probe walkers run between chunks is a
relaxed atomic read, so an idle collision queue costs no lock traffic.

### 3.3 Per-thread wind-down state

Each worker thread carries a `state` that the comm thread *writes to ask* and the
thread itself *writes to answer*; a thread always announces its own completion.

```
walker:    RUNNING --[comm]--> HOLD --[self]--> HELD --[comm]--> DRAIN --[self]--> QUIESCENT
inserter:  RUNNING -----------------------------------[comm]--> DRAIN --[self]--> QUIESCENT
```

| state | walker | inserter |
|---|---|---|
| `RUNNING` | walks chunks, ships DPs, retires collision candidates between chunks | probes everything in its queue |
| `HOLD` | "stop producing points": acknowledges by writing `HELD` | (not used) |
| `HELD` | no new points; keeps retiring collision candidates | (not used) |
| `DRAIN` | "no candidate will ever be added": empties `coll_q`, publishes `n_dp`, writes `QUIESCENT`, returns | "no point will ever be delivered": empties its queue, writes `QUIESCENT`, returns |
| `QUIESCENT` | done for this round | done for this round |

A walker reads its state **once per chunk** (`chunk_size` iterations of `vmixf`,
default 64, plus whatever candidates it retires before the chunk), which is why
`HELD` exists: until a walker has acknowledged, an empty walker queue is not proof
that it is idle.  An inserter checks its state **before** re-testing emptiness, so a
point pushed between the two checks is still probed before it goes quiet.  All state
loads are acquire, all stores release.  Every thread resets itself to `RUNNING` at the
top of the next round.

### 3.4 Statistics and the golden pair

- `n_dp` (walkers), `n_probe`, `n_drop_coll` (inserters), `n_drop_walkerq` (walkers)
  are running totals for the round, stored relaxed by the owner (a walker publishes
  `n_dp` after every chunk) and read relaxed by the comm thread when it builds a
  progress report.  The comm thread zeroes them in the epilogue.
- `Counters ctr` is written by its thread with **no synchronisation at all** and is
  read by the comm thread only after the end-of-round OpenMP barrier.  It is merged
  across threads, reduced across nodes, and reset by the comm thread (§4.6).
- `RoundState::set_golden(i, x0, x1)` takes a mutex, keeps the first triple, and
  raises `golden_found` with a release store.  The comm thread polls the flag (acquire)
  every turn of its loop during the round and forwards it as a golden report.
  Neither `golden_found` nor `golden_sent` is ever reset: a golden pair ends the search.

## 4. Lifecycle and synchronization steps

### 4.1 Startup (before any round)

1. The driver calls `MPI_Init_thread(MPI_THREAD_FUNNELED)`; `run_engine()` refuses a
   lower level.  Every rank must use the same PRNG seed (the driver broadcasts it).
2. `Parameters` is built on every rank from identical inputs (rank and size come from
   `MPI_Comm_rank/size`).
3. The test-vector `MPI_Bcast` (§2.3) asserts every rank iterates the same function.
4. `PcsNode` is constructed: `InBuffers` posts its `n_in_buffers` receives,
   `ControlChannel` attaches the `MPI_Bsend` buffer and posts the control receive.
   Nothing may `MPI_Bsend` before this, and nothing does.
5. The OpenMP team starts.  Each thread pins itself to `thread_cpu[tid]`.
   **OpenMP barrier**, then the master prints the thread-to-CPU map on rank 0.

### 4.2 Round start

1. Every thread sets its own `state = RUNNING`.
2. **Master only** (thread 0): rank 0 draws `i` and `root_seed` and reads
   `controller.stop`; **`MPI_Bcast` of the 3-word round header** from rank 0.  Every
   rank stores it into `RoundState`.  On rank 0, `controller.begin_round()` resets the
   per-round tallies and `n_active_nodes = n_nodes`, unless `stop` is set.
3. **OpenMP barrier**: publishes `RoundState` to every thread.
4. If `stop`, every thread leaves the round loop (§4.7).  Otherwise thread 0 enters
   `comm_round()`, inserters `inserter_thread()`, walkers `walker_thread()`.

### 4.3 Steady state

**Comm thread**, one turn of its loop, forever until `NEW_VERSION` arrives:

1. `route_walker_queues()`: pop up to 64 DPs from each walker queue; local points go
   to `deliver_local()`, remote ones to `outbuf.push()`.
2. `outbuf.poll()`: `MPI_Testsome` on the outgoing sends, releasing finished buffers.
3. `poll_incoming()`: `MPI_Testsome` on the receive pool; scatter completed buffers to
   the inserter queues, count sentinels, repost.
4. `service_control()`: `MPI_Test` the control receive.  An assignment clears
   `awaiting_assignment` and, if it is `NEW_VERSION`, sets `got_new_version`.  A
   report (rank 0 only) goes to `controller.handle_report()`, whose reply, if any, is
   `MPI_Bsend`-ed back to the report's source.
5. If `golden_found` and not yet sent: `MPI_Bsend` a golden report to rank 0.
6. Every 256 turns, if no report is outstanding and the pacing rule (§2.2) says so:
   read the walkers' and inserters' published tallies, `MPI_Bsend` a progress report
   with the deltas, set `awaiting_assignment`.
7. If `got_new_version`: leave the loop and run the drain (§4.5).

Nothing in this loop blocks.

**Walker**, one turn: read `state`; retire up to `coll_per_chunk` collision candidates
(0 = all); walk one chunk, pushing each DP found into its SPSC queue (drop if full)
and restarting the chain; add the chunk's evaluations to `ctr`; publish `n_dp`.

**Inserter**, one turn: pop up to 64 DPs, probe each into its shard
(`PcsDict::pop_insert`), push a `CollisionCandidate` for every hit (drop if full);
then read `state`; if nothing to do, `cpu_relax()`.

**Controller** (inside `handle_report` on rank 0): a golden report records the
solution (first one wins) and raises `stop`; a progress report adds its deltas, and is
answered `NEW_VERSION` if `stop` is raised or `ndp >= points_per_version`, else
`KEEP_GOING`.  The live one-line display is refreshed at most every 0.5 s.

### 4.4 Point flow, end to end

```
walker finds a DP
   |  SPSC
   v
comm thread (sender):  node = x % n_nodes
   |                                    \
   | local                               \  remote: OutBuffers, MPI_Isend TAG_POINTS
   v                                      v
deliver_local()   <-----------------  comm thread (receiver): poll_incoming()
   |  slot = (x / n_nodes) % inserters_per_node, SPSC
   v
inserter:  key = x / n_inserters, probe the shard
   |  hit: CollisionCandidate, mutex
   v
coll_q  -->  walker: walk both trails, locate the collision, test the pair
                |  golden
                v
          RoundState::set_golden  -->  comm thread  -->  golden report to rank 0
```

### 4.5 End of round: the drain

Entered by the comm thread when its node has been told `NEW_VERSION`.  The order is
what makes the round airtight: no thread declares itself finished while something
can still arrive for it.

| step | comm thread does | guarantee once it passes |
|---|---|---|
| 1 | `state(WALKER) := HOLD` | walkers will stop producing at their next chunk boundary; they keep resolving collisions |
| 2 | loop { check every walker is `HELD`; `route_walker_queues()`; `outbuf.poll()`; `poll_incoming()`; `service_control()` } until all walkers are `HELD`, the last pass moved nothing, and every walker queue has `head == tail` | every DP this node produced this round has been routed (delivered locally or handed to `OutBuffers`) |
| 3 | loop { `outbuf.flush_poll()`; `poll_incoming()`; `service_control()` } until nothing is left to send and every `Isend` has completed | every DP bound for another node has left this node.  Receiving continues throughout, so the peers we wait on can complete their sends too |
| 4 | `outbuf.send_sentinels()` | every node (self included) will learn that we are done sending; the sentinel cannot overtake the data it follows |
| 5 | loop { `poll_incoming()`; `service_control()` } until `n_sentinels == n_nodes` | every node has finished sending to us, and everything they sent has been scattered to the inserter queues.  Rank 0 keeps answering reports here: the other nodes reach *their* step 4 only after it tells them `NEW_VERSION` |
| 5b | rank 0 only: loop { `service_control()` } until `controller.n_active_nodes == 0` | rank 0 does not leave the round while a node is still waiting for its assignment (which would strand it in front of the collectives).  Already implied by step 5 in practice, since a node only sends its sentinel after `NEW_VERSION` |
| 6 | `state(INSERTER) := DRAIN`; loop { `service_control()` } until all inserters are `QUIESCENT` | every DP of the round has been probed; no new collision candidate can appear |
| 7 | `state(WALKER) := DRAIN`; loop { `service_control()` } until all walkers are `QUIESCENT` | the collision queue is empty; every candidate of the round has been resolved.  Doing this last is what keeps a golden pair found on the very last candidate from being lost |

The control channel is serviced in every waiting loop because rank 0 hands itself
`NEW_VERSION` like everyone else and drains like everyone else; if it stopped
answering during its own drain, the other nodes would never be told to stop.

### 4.6 Epilogue

1. **OpenMP barrier**: every worker of this node has returned from its round function
   and all their `Counters` are now safe to read.
2. Thread 0: merge the `n_threads` `Counters` into one, pack the `ST_*` words,
   **`MPI_Reduce` (SUM) to rank 0**, then **`MPI_Reduce` (MAX) of the HyperLogLog
   registers**.  On rank 0, `controller.end_round()` prints the round report, folds
   the totals, bumps `nround`, and raises `stop` if `nround >= max_versions`.
   Then thread 0 resets every thread's `ctr` and published tallies, `outbuf.n_dropped`,
   `n_drop_inserterq` and `inbuf.n_sentinels`.
   Meanwhile, each inserter thread zeroes its own dictionary shard (`PcsDict::flush`).
3. **OpenMP barrier**: the next round's `Bcast` may not start before every shard is
   clean and the tallies are cleared.

### 4.7 Termination

When the round header carries `stop`, every thread of every rank leaves the loop
after the barrier of §4.2 step 3 (no round is armed).  Thread 0 of rank 0 packs the
solution (`found, i, x0, x1`), prints the final line, and **`MPI_Bcast`s the 4-word
answer**; every rank then cancels its posted receives (`InBuffers::shutdown`,
`ControlChannel::shutdown`, which also detaches the `Bsend` buffer) and returns from
`run()`.  The driver calls `MPI_Finalize`.

## 5. Guarantees and loss semantics

**What may be lost.**  Distinguished points and collision candidates are
loss-tolerant: losing one costs the work that produced it and nothing else, so every
DP/candidate channel drops on overflow rather than blocking, and every drop is
tallied and shown in the round report (`DROPPED ... walker-queue / output-buffer /
inserter-queue / collision-queue`).

| channel | may drop | counted in |
|---|---|---|
| walker -> comm SPSC | yes | `n_drop_walkerq` |
| comm -> remote node (`OutBuffers`) | yes | `n_dropped` (`REP_DROP_OUT`) |
| comm -> inserter SPSC | yes | `n_drop_inserterq` |
| inserter -> walkers (`coll_q`) | yes | `n_drop_coll` |
| control channel (reports, assignments) | **never** | buffered send, buffer sized for the worst case |
| sentinels | **never** | same buffer |
| collectives | n/a | |

**What is enforced per round.**

- Every DP produced in round `r` is either dropped (and counted) or probed into the
  round-`r` dictionary, on the shard that owns it, before that shard's inserter goes
  quiescent (drain steps 2 through 6).  No round-`r` point is ever probed in round
  `r+1`: nothing is delivered after the last sentinel, and the queues are empty when
  the inserters stop.
- Every collision candidate of round `r` is resolved by a round-`r` walker before the
  walkers go quiescent (drain step 7); `service_collision` asserts `c.i == round.i`.
- The dictionary is empty at the start of every round (inserters flush between the
  two epilogue barriers).
- A node reports at most once between two assignments, so the controller counts each
  node once when it switches version (see §6 for the one window).
- The controller never initiates a message; it only answers, so a round can only end
  through the report/assignment handshake.

**Why it does not deadlock.**

- The comm thread never blocks on point-to-point MPI: every send is `Isend` (tested,
  never waited) or `Bsend` (completes locally); every receive is an always-posted
  `Irecv` that is only ever `Test`ed.
- While waiting for its own sends to complete (step 3) or for the peers' sentinels
  (step 5), a node keeps receiving and keeps servicing the control channel, so two
  nodes waiting on each other both make progress.
- Rank 0 services reports during every phase of its own drain, so it can always hand
  the remaining nodes their `NEW_VERSION`.
- The only blocking MPI calls (the round-start `Bcast` and the epilogue reductions)
  are reached by a rank only after its drain, i.e. after every node has been told
  `NEW_VERSION` and has confirmed with a sentinel; so every rank reaches them.
- Inside a node, no thread ever waits on a queue; the only waits are the comm thread
  spinning on `state` values that the worker threads set for themselves at chunk
  boundaries, and the three OpenMP barriers per round, which every thread reaches
  unconditionally once its round function returns.

## 6. Known windows

Two places where the protocol is slightly weaker than the rules above suggest.  Both
are benign for the search itself; they are listed so nobody rediscovers them.

- **A golden pair found during the drain is reported one round late.**  The comm
  thread forwards `golden_found` only inside its steady-state loop (§4.3 step 5);
  the drain loops service the control channel but do not check the flag.  Walkers
  keep resolving candidates through drain steps 1 through 7, so a pair found there
  is sent at the first turn of the *next* round's loop, after which the controller
  raises `stop`, hands out `NEW_VERSION` to everyone, and the search ends one round
  later than it could.  If that drain belonged to the last permitted round
  (`max_versions`), no next round starts and the pair is never reported.
- **One extra progress report can slip out after `NEW_VERSION`.**  Within one turn of
  the comm loop, `service_control()` (which clears `awaiting_assignment` and sets
  `got_new_version`) runs *before* the periodic-report check, and the loop only
  breaks at the end of the turn.  On the 1-in-256 turn where the check fires and the
  pacing rule allows it, a report leaves after the node has already been counted
  out.  The controller then either decrements `n_active_nodes` a second time (harmless:
  step 5, not 5b, is the real gate) or, if the report is polled after the next
  `begin_round()`, credits a few stale DPs to the new round and answers with a
  `KEEP_GOING` that clears the node's `awaiting_assignment` flag one report early.
