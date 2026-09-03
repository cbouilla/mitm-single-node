# Communication protocol

Everything below is implemented in `include/engine.hpp` (the comm thread and the
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
and **only thread 0 ever calls MPI** (its code is the `CommThread` class); the
workers never see an `MPI_` symbol.
Worker threads communicate with thread 0 through per-thread
Single-Producer-Single-Consumer (SPSC) queues, one `state` atomic each and their
plain-`u64` tallies, and with each other only through the collision queue.

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
**dropped** (the comm thread tallies it, `DROP_OUT`), never blocked on.  A buffer is
never sent empty.

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

### 2.2 The control channel: `TAG_END_ROUND`, `TAG_REPORT`, `TAG_SOLUTION`

Three tags, one per message kind, and one always-posted `MPI_Irecv` per tag on the
rank that consumes it:

| tag | from | to | words | posted by |
|---|---|---|---|---|
| `TAG_END_ROUND` | rank 0 | every node | 0 | `ControlChannel`, on every rank, with source `0` (no wildcard) |
| `TAG_REPORT` | every node | rank 0 | `N_COUNTERS` | `Controller`, rank 0 only, `MPI_ANY_SOURCE` |
| `TAG_SOLUTION` | every node | rank 0 | `SOL_NWORDS` | `Controller`, rank 0 only, `MPI_ANY_SOURCE` |

Which receive completed says what the message is, so each buffer is sized exactly for
its message and nothing is ever inferred from a length.  There is no `MPI_ANY_TAG`
anywhere, which is what keeps these receives from ever matching `TAG_POINTS` traffic
bound for rank 0's node role.  The report and solution receives are serviced by the
controller itself (`Controller::service()`); the node's comm code only ever sees the
end-of-round signal.  Rank 0 talks to itself through MPI like any other node.  Every
send is `MPI_Bsend`: it completes locally, so the comm thread never blocks and there is
no request to track.  The attached buffer is per process and holds
`2 * n_nodes + bsend_slack` report-sized slots.  The bounded part is an end-of-round
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

**Progress report** (node -> rank 0, `TAG_REPORT`, `N_COUNTERS == 16` words).  The
node's `Counters` array -- every thread's tallies summed (§3.4), in `enum counter`
order (`counters.hpp`) -- as **deltas since the node's previous report**; the
controller sums them into `reported[]`.  Only two of them act during the round:
`N_DP` closes it and `N_PROBE` feeds the live display.  The rest ride along because
one layout serves the per-thread tallies, the report and the end-of-round reduction
alike.

| index | meaning | tallied by |
|---|---|---|
| `N_EVAL` | evaluations of the mixing function, walking or resolving | walkers |
| `N_DP` | distinguished points found | walkers |
| `N_POINTS_TRAILS` | sum of the lengths of the trails that reached a DP | walkers |
| `N_COLLISIONS` | collisions located | walkers |
| `COLLIDING_LEN_MIN` | sum of the shorter trail length of each colliding pair | walkers |
| `COLLIDING_LEN_MAX` | ... and of the longer one | walkers |
| `BAD_DP` | trail gave up before reaching a distinguished point | walkers |
| `BAD_COLLISION` | the two trails "collide" on the same value | walkers |
| `BAD_WALK_ROBINHOOD` | one trail is a suffix of the other | walkers |
| `BAD_WALK_NONCOLLIDING` | dictionary false positive: the trails never meet | walkers |
| `DROP_WALKERQ` | points dropped: walker queue to comm full | walkers |
| `N_PROBE` | dictionary probes retired | inserters |
| `BAD_PROBE` | probe missed: slot empty, or held a different key | inserters |
| `DROP_COLL` | candidates dropped: collision queue full | inserters |
| `DROP_OUT` | points dropped: outgoing MPI buffer busy | comm thread |
| `DROP_INSERTERQ` | points dropped: inserter queue full | comm thread |

Pacing: the comm thread considers reporting every 256 turns of its poll loop, and
sends a report if either `ping_delay` seconds (default 0.1) have elapsed since the
last one, or the node has found `points_per_version / (reports_per_round * n_nodes)`
new DPs since the last one (so a short round cannot overshoot `beta * w`).  A report
is **one-way**: nothing is answered, and a node keeps reporting until its end-of-round
signal arrives.

**Solution** (node -> rank 0, `TAG_SOLUTION`, `SOL_NWORDS == 3` words).  `SOL_I`,
`SOL_X0`, `SOL_X1`: the mixing-function index and the two colliding points, exactly
the layout of `RoundState::golden`.  Sent once per node (`golden_sent`) as soon as the
comm thread sees `RoundState::golden_found`.  It is **one-way**.  The controller keeps
the first one and raises `stop`.

**End of round** (rank 0 -> every node, `TAG_END_ROUND`, **zero-length**).  Meaning:
"stop producing points and drain: the round is over".  Sent by
`Controller::close_round()` to every node, rank 0 included, as soon as `service()`
sees `stop` raised or the summed `N_DP` of the round's reports at or above
`points_per_version`.  It is the only message the controller ever sends, it is
unsolicited, and it is sent **exactly once per round** (`round_closed`): the receive
that matches it names source 0 and this tag, so successive signals are non-overtaking
and a node consumes exactly one per round -- it is its only way out of the steady
state -- and a second one would be consumed in the *next* round and end it at once.
The controller only knows what has been reported, so a node that has not reported
cannot push `ndp` over the threshold; that is why reports are also paced by DP count.

### 2.3 Collectives

All on `mpi_comm`, all issued by thread 0, all `MPI_UINT64_T` unless stated.

| when | call | payload | root |
|---|---|---|---|
| driver `mitm::init`, if the seed was not given | `MPI_Bcast` | 1 word: PRNG seed (this is in `examples/driver.hpp`, not the library) | 0 |
| `run()`, before the team exists | `MPI_Bcast` | 3 words: `x`, `y`, `mixf(x, y)`; every rank asserts it computes the same value | 0 |
| start of every round | `MPI_Bcast` | 3 words: `i` (function version), `root_seed`, `stop` | 0 |
| end of every round | `MPI_Reduce`, `MPI_SUM` | `N_COUNTERS == 16` words: the round's merged `Counters` array, in `enum counter` order (the layout of a report, totals instead of deltas) | 0 |
| end of every round | `MPI_Reduce`, `MPI_MAX` | 65 536 `MPI_UINT8_T`: HyperLogLog registers of this round's collisions | 0 |
| after the last round | `MPI_Bcast` | 4 words: `found`, `i`, `x0`, `x1` | 0 |

These collectives are the **only blocking MPI calls** the engine makes after
startup.  The per-round ones are reached only after every node has finished its drain
(§4.5), so no point-to-point traffic of the round can still be pending when they run.

## 3. Channels inside a node

| channel | producer | consumer | type | capacity (`Options`) | when full |
|---|---|---|---|---|---|
| `ctx[R+1+w]->q` | walker `w` | comm | `SPSCQueue` of `DP` | `walker_queue_capacity` (1024) | walker drops the DP, tallies `DROP_WALKERQ` |
| `ctx[1+r]->q` | comm | inserter `r` | `SPSCQueue` of `DP` | `inserter_queue_capacity` (4096) | comm drops the DP, tallies `DROP_INSERTERQ` |
| `coll_q` | every inserter | every walker | `CollisionQueue` (mutex, bounded) | `coll_queue_capacity` (8192) | inserter drops the candidate, tallies `DROP_COLL` |
| `ThreadContext::state` | comm and the thread itself | the other side | `atomic<int>` | | see §3.3 |
| `ThreadContext::ctr` | the thread | comm | plain `Counters`: `u64[N_COUNTERS]` + HLL registers, no atomics, see §3.4 | | |
| `RoundState::i`, `root_seed`, `stop` | comm (thread 0) | everyone | plain `u64`, published by an OpenMP barrier | | |
| `RoundState::golden` | any walker | comm | mutex + `atomic<bool>` flag | | first one wins |

Each SPSC queue belongs to the `ThreadContext` of its worker side (the walker, or
the inserter), which builds it once that thread is pinned (§4.1).  SPSC capacities are rounded up to a
power of two.  Every queue is **non-blocking on both sides**: nobody ever waits for
room or for data, they poll and move on.

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
| `DRAIN` | "no candidate will ever be added": empties `coll_q`, writes `QUIESCENT`, returns | "no point will ever be delivered": empties its queue, writes `QUIESCENT`, returns |
| `QUIESCENT` | done for this round | done for this round |

A walker reads its state **once per chunk** (`chunk_size` iterations of `vmixf`,
default 64, plus whatever candidates it retires before the chunk), which is why
`HELD` exists: until a walker has acknowledged, an empty walker queue is not proof
that it is idle.  An inserter checks its state **before** re-testing emptiness, so a
point pushed between the two checks is still probed before it goes quiet.  All state
loads are acquire, all stores release.  Every thread resets itself to `RUNNING` at the
top of the next round.

### 3.4 Statistics and the golden pair

- Every tally of a thread lives in its `Counters ctr`: one `u64[N_COUNTERS]` indexed
  by `enum counter`, plus the HyperLogLog registers.  The owner writes it with **no
  synchronisation at all** -- plain increments, no atomics -- and the comm thread's
  own drops go into `ctx[0]->ctr` the same way.  The comm thread reads them twice.
  During the round, `CommThread::snapshot()` sums the arrays of every thread after a
  `#pragma omp flush`, to build a progress report: it sees whatever has reached
  memory, so a count may lag by a unit or two, which is fine for a report.  After the
  end-of-round OpenMP barrier the values are exact: they are merged across threads,
  reduced across nodes, and reset by the comm thread (§4.6).
- `RoundState::set_golden(i, x0, x1)` takes a mutex, keeps the first triple, and
  raises `golden_found` with a release store.  The comm thread polls the flag (acquire)
  every turn of its loop during the round and forwards it as a solution message.
  Neither `golden_found` nor `golden_sent` is ever reset: a golden pair ends the search.

## 4. Lifecycle and synchronization steps

### 4.1 Startup (before any round)

1. The driver calls `MPI_Init_thread(MPI_THREAD_FUNNELED)`; `run()` refuses a lower
   level.  Every rank must use the same PRNG seed (the driver broadcasts it).
2. `run()`, on the main thread: `Parameters` is built on every rank from identical
   inputs (rank and size come from `MPI_Comm_rank/size`); the test-vector `MPI_Bcast`
   (§2.3) asserts every rank iterates the same function; rank 0 prints the banner.
   Then what the threads share is built, as locals of `run()`: the per-thread table
   `ctx` (sized, not filled), the collision queue, the `RoundState`, and the
   `CommThread` -- thread 0's object, and the main thread *is* thread 0 of the team
   to come.  Its constructor posts every receive of the engine: `InBuffers` posts its
   `n_in_buffers` receives, `ControlChannel` attaches the `MPI_Bsend` buffer and posts
   the end-of-round receive, and on rank 0 `Controller` posts the report and solution
   receives.  Nothing may `MPI_Bsend` before this, and nothing does.  Nothing big is
   allocated before the team exists: the MPI buffers are a few tens of kB.
3. The OpenMP team starts.  Each thread pins itself to `thread_cpu[tid]`, **then**
   builds its `ThreadContext`, which holds what it owns: for a walker its SPSC queue,
   for an inserter its SPSC queue and its dictionary shard (`PcsDict`, `w_shard`
   slots, zero-filled).  Allocating after pinning is deliberate: under Linux's
   first-touch policy a page lands on the NUMA node of the CPU that first writes it,
   so this puts every object on its owner's node.  **OpenMP barrier** publishes the
   slots; the comm thread reads them from its first round on.

### 4.2 Round start

1. Every thread sets its own `state = RUNNING`.
2. **Thread 0 only**, `CommThread::begin_round()`: rank 0 draws `i` and `root_seed`
   and reads `controller.stop`; **`MPI_Bcast` of the 3-word round header** from
   rank 0.  Every rank stores it into `RoundState`.  On rank 0,
   `controller.begin_round()` resets the per-round tallies and `round_closed`, unless
   `stop` is set.
3. **OpenMP barrier**: publishes `RoundState` to every thread.
4. If `stop`, every thread leaves the round loop (§4.7).  Otherwise thread 0 enters
   `CommThread::comm_round()`, inserters `inserter_thread()`, walkers `walker_thread()`.

### 4.3 Steady state

**Comm thread**, one turn of its loop, forever until the end-of-round signal arrives:

1. `route_walker_queues()`: pop up to 64 DPs from each walker queue; local points go
   to `deliver_local()`, remote ones to `outbuf.push()`.
2. `outbuf.poll()`: `MPI_Testsome` on the outgoing sends, releasing finished buffers.
3. `poll_incoming()`: `MPI_Testsome` on the receive pool; scatter completed buffers to
   the inserter queues, count sentinels, repost.
4. `service_control()`: `MPI_Test` the end-of-round receive; if it completed, set
   `round_over` (and repost it for the next round).  On rank 0, then
   `controller.service()`: drain the solution receive (each one goes to
   `handle_solution()`) and the report receive (each one goes to `handle_report()`),
   then close the round (§2.2) if it is not closed yet and `stop` is raised or
   `ndp >= points_per_version`.
5. If `golden_found` and not yet sent: `MPI_Bsend` the solution to rank 0.
6. Every 256 turns: `snapshot()` the node's tallies (§3.4) and, if the pacing rule
   (§2.2) says so, `MPI_Bsend` a progress report with the deltas since the previous one.
7. If `round_over`: leave the loop and run the drain (§4.5).

Nothing in this loop blocks.

**Walker**, one turn: read `state`; retire up to `coll_per_chunk` collision candidates
(0 = all); walk one chunk, pushing each DP found into its SPSC queue (drop if full)
and restarting the chain, tallying into `ctr` as it goes; add the chunk's evaluations
to `ctr`.

**Inserter**, one turn: pop up to 64 DPs, probe each into its shard
(`PcsDict::pop_insert`), push a `CollisionCandidate` for every hit (drop if full);
then read `state`; if nothing to do, `cpu_relax()`.

**Controller** (inside `service()` on rank 0): a solution is recorded (first one wins)
and raises `stop`; a progress report adds its deltas; then, once per round, the round
is closed -- `TAG_END_ROUND` to every node -- as soon as `stop` is raised or
`ndp >= points_per_version`.  The live one-line display is refreshed at most every
0.5 s.

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
          RoundState::set_golden  -->  comm thread  -->  solution to rank 0
```

### 4.5 End of round: the drain

Entered by the comm thread when its end-of-round signal has arrived.  The order is
what makes the round airtight: no thread declares itself finished while something
can still arrive for it.

| step | comm thread does | guarantee once it passes |
|---|---|---|
| 1 | `state(WALKER) := HOLD` | walkers will stop producing at their next chunk boundary; they keep resolving collisions |
| 2 | loop { check every walker is `HELD`; `route_walker_queues()`; `outbuf.poll()`; `poll_incoming()`; `service_control()` } until all walkers are `HELD`, the last pass moved nothing, and every walker queue has `head == tail` | every DP this node produced this round has been routed (delivered locally or handed to `OutBuffers`) |
| 3 | loop { `outbuf.flush_poll()`; `poll_incoming()`; `service_control()` } until nothing is left to send and every `Isend` has completed | every DP bound for another node has left this node.  Receiving continues throughout, so the peers we wait on can complete their sends too |
| 4 | `outbuf.send_sentinels()` | every node (self included) will learn that we are done sending; the sentinel cannot overtake the data it follows |
| 5 | loop { `poll_incoming()`; `service_control()` } until `n_sentinels == n_nodes` | every node has finished sending to us, and everything they sent has been scattered to the inserter queues.  Every node was sent its end-of-round signal before rank 0 could even enter its drain, so nobody is waiting on rank 0 here; it keeps digesting reports so the round's tallies cover what the others produced meanwhile |
| 6 | `state(INSERTER) := DRAIN`; loop { `service_control()` } until all inserters are `QUIESCENT` | every DP of the round has been probed; no new collision candidate can appear |
| 7 | `state(WALKER) := DRAIN`; loop { `service_control()` } until all walkers are `QUIESCENT` | the collision queue is empty; every candidate of the round has been resolved.  Doing this last is what keeps a golden pair found on the very last candidate from being lost |

The control channel is serviced in every waiting loop: on rank 0 so that the reports
and solutions of nodes still in their steady state are digested (the round's tallies,
and `stop` for the next round header), and on every rank because in steps 6 and 7 it is
the only MPI call, whose progress our sentinels and rank 0's signals need to actually
leave.  No node's liveness depends on it: every signal of the round has been sent
before rank 0 can enter its own drain.

### 4.6 Epilogue

1. **OpenMP barrier**: every worker of this node has returned from its round function
   and all their `Counters` are now exact and safe to read.
2. Thread 0, `CommThread::end_round()`: merge the `n_threads` `Counters` into one,
   **`MPI_Reduce` (SUM) of its `N_COUNTERS` words to rank 0**, then **`MPI_Reduce`
   (MAX) of its HyperLogLog registers**.  On rank 0, `controller.end_round()` prints
   the round report from the reduced values, folds them into its all-time `Counters`,
   bumps `nround`, and raises `stop` if `nround >= max_versions`.  Then thread 0 resets
   every thread's `ctr` and `inbuf.n_sentinels`.  Meanwhile, each inserter thread
   zeroes its own dictionary shard (`PcsDict::flush`).
3. **OpenMP barrier**: the next round's `Bcast` may not start before every shard is
   clean and the tallies are cleared.

### 4.7 Termination

When the round header carries `stop`, every thread of every rank leaves the loop
after the barrier of §4.2 step 3 (no round is armed).  Thread 0 then runs
`CommThread::finish()`: on rank 0 it packs the solution (`found, i, x0, x1`) and
prints the final line; every rank **`MPI_Bcast`s the 4-word answer**, then cancels
its posted receives (`InBuffers::shutdown`, `Controller::shutdown`, a no-op off
rank 0, then `ControlChannel::shutdown`, which also detaches the `Bsend` buffer).
`run()` returns the answer on every rank.  The driver calls `MPI_Finalize`.

## 5. Guarantees and loss semantics

**What may be lost.**  Distinguished points and collision candidates are
loss-tolerant: losing one costs the work that produced it and nothing else, so every
DP/candidate channel drops on overflow rather than blocking, and every drop is
tallied and shown in the round report (`DROPPED ... walker-queue / output-buffer /
inserter-queue / collision-queue`).

| channel | may drop | counted in |
|---|---|---|
| walker -> comm SPSC | yes | `DROP_WALKERQ` |
| comm -> remote node (`OutBuffers`) | yes | `DROP_OUT` |
| comm -> inserter SPSC | yes | `DROP_INSERTERQ` |
| inserter -> walkers (`coll_q`) | yes | `DROP_COLL` |
| control channel (end-of-round signals, reports, solutions) | **never** | buffered send, buffer sized for the bounded traffic plus slack (§2.2); a full buffer is a fatal error, not a drop |
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
- The controller sends exactly one message per round, the end-of-round signal to every
  node, and answers nothing; a round can end no other way, and a node leaves the
  steady state exactly once per round.

**Why it does not deadlock.**

- The comm thread never blocks on point-to-point MPI: every send is `Isend` (tested,
  never waited) or `Bsend` (completes locally); every receive is an always-posted
  `Irecv` that is only ever `Test`ed.
- While waiting for its own sends to complete (step 3) or for the peers' sentinels
  (step 5), a node keeps receiving and keeps servicing the control channel, so two
  nodes waiting on each other both make progress.
- Rank 0 closes the round for every node in one go, and observes its own signal only on
  a later turn, so by the time it enters its own drain every node has been sent its
  signal.
- The only blocking MPI calls (the round-start `Bcast` and the epilogue reductions)
  are reached by a rank only after its drain, i.e. after every node has received the
  end-of-round signal and has confirmed with a sentinel; so every rank reaches them.
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
  raises `stop`, closes the round for everyone, and the search ends one round
  later than it could.  If that drain belonged to the last permitted round
  (`max_versions`), no next round starts and the pair is never reported.
- **Reports are one-way, so the round's `reported[]` tallies are approximate.**  A
  node keeps reporting until its end-of-round signal is matched, and the points it
  finds between its last report and that moment are never reported; on top of that a
  snapshot reads the workers' tallies without synchronisation (§3.4).  Only the live
  display and the decision to close the round rest on them; the reductions of §4.6
  count everything exactly, and the round report prints those.  Further, a report and
  the same node's later sentinel match different receives, so the MPI standard does
  not order them: a report could in principle be matched only after rank 0 has run
  the epilogue and reset `reported[]`, and be credited to the next round, which could
  then close early by at most one report per node.  Open MPI matches messages per
  source in sequence order and `service()` drains the report receive to empty on
  every call, so this does not happen in practice.
