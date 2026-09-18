# Communication protocol

**Two engines live in this tree, and this document covers both.**  Both are built on the
Router (`include/router/`, its interface specified in `router.3`), both run one MPI rank per
node and one OpenMP team per rank, and neither shares a line of transport with the other.

- **§1, the direct engine** (`include/direct/`): the exhaustive meet-in-the-middle.  Its
  round is deterministic, so it needs no control channel at all.
- **§2, PCS** (`include/pcs/`): van Oorschot-Wiener parallel collision search over a
  dictionary of trail endpoints.  Its round ends on a *global* count, which no node can tell
  by itself, so it keeps a controller on rank 0 and a control channel to reach it.

All either engine asks of the communicator is that nothing else on it carries the tags it
owns.  Names in code font are the ones used in the code.

## 1. The direct engine, on the Router

`mitm::direct` (`include/direct/`): the exhaustive meet-in-the-middle, the baseline PCS will be
measured against.  `w' = fill * w` entries per round (`--fill`, default 0.5),
`R = ceil(2^n / w')` rounds; `--nrounds` caps `R`, and the search is then not exhaustive -- the
banner and the last line say so.  Nothing here is shared with §2.

### 1.1 The team

One MPI rank per node, `MPI_THREAD_FUNNELED`, one OpenMP team of `1 + I + W` threads per rank
(`I = --dicts-per-node`, `W = --producers-per-node`, 0 == fill the affinity mask).  Every thread
calls `Router_Init(role, ROUTER_GROUP_AUTO, comm, ROUTER_TAG, &opts.router)` once,
inside the one parallel region, and keeps its handle for the team's whole life.

| thread | Router role | does |
|---|---|---|
| 0 | service | all the MPI: `Router_Progress` until `Router_Test_quiescent`, then `Router_Stats`, the epilogue, the printing |
| `1..I` | receiver | owns one dictionary shard: inserts or probes what the Router delivers to it |
| `I+1..I+W` | sender | evaluates the phase's function on its piece of the domain and pushes every image |

**The Router owns placement** (autopilot): it pins every thread, forms its groups and reports the
plan and the measured layout itself.  The engine passes `ROUTER_GROUP_AUTO`, never pins a thread and
never queries hwloc; `--no-bind` is `opts.router.pin = false`, needed when two ranks share a host.
The engine ignores the groups: a match is resolved where it is found, so no producer needs to be
paired with a dict thread.  What it does ask the Router is **the team's shape**, through
`Router_size(rt, role, scope)` and `Router_rank(rt, scope)`: the shard count for the routing
(§1.3), the producer count and a producer's global rank for its piece of the span, the node count
for the epilogue's record array (§1.4), the node's dict threads for the shard report.
`direct::Params` holds the dictionary's slots, the domain and the round count, nothing about the
team.  A dict thread builds its shard right after `Router_Init`, i.e. once
pinned, so the zero-fill is the first touch of every page.  That shard is one `mmap`, asked for on
2 MB pages (`MAP_HUGETLB`) and taken on ordinary ones when the kernel has no huge page reserved; the
length is rounded up to a whole huge page, so a run holds up to 2 MB per shard more than `--ram`
allows, and what the kernel really granted is read back from `/proc/self/smaps` and reported.

### 1.2 One direct round is two Router rounds

The phase sequence is `(0, FILL)`, `(0, PROBE)`, `(1, FILL)`, ... and it is **deterministic**: every
thread of every node computes it from a local counter, so there is no round header, no broadcast and
no control channel.

| phase | a producer evaluates | on | a dict thread | at the end of the phase |
|---|---|---|---|---|
| `FILL` | `f` | chunk `round` of the domain: `[round * w', min((round + 1) * w', 2^n))` | inserts `(key, val)` | keeps the shard |
| `PROBE` | `g` (`f` for a collision problem) | the whole domain | probes `key`; resolves every candidate on the spot | `clear()`s the shard |

**The phase boundary is `Router_Reset`**, and that is the barrier the algorithm needs: it returns on
a thread only after every node has left the round, so every entry of a `FILL` phase is in its shard
before the first probe of the `PROBE` phase is pushed.  Its two team barriers and its `MPI_Barrier`
are inside it; the engine writes none of them.

One phase, on each thread:

1. zero its own tallies;
2. its role's work, above -- a producer ends with `Router_Close`, a dict thread loops until
   `Router_Pop` finds nothing and `Router_Test_drained` is true, thread 0 turns the Router until the
   node is quiescent and then reads `Router_Stats`, which `Router_Reset` would clear;
3. `Router_Reset`, which every thread of every node calls;
4. thread 0 alone: the epilogue (§1.4);
5. `#pragma omp barrier`, the engine's own: it publishes the epilogue's verdict and lets the next
   phase zero the tallies thread 0 has just read;
6. every thread leaves the loop if `found` is set or the domain is exhausted -- a test it makes for
   itself -- or starts the next phase.

### 1.3 Points

A point is the two words the Router carries: `key` is `murmur64(image)`, `val` the preimage.

**Routing.**  `dest = (key * n_dicts) >> 64`, a multiply-shift on the key's top bits, where
`n_dicts` is `Router_size(rt, ROUTER_RECEIVER, ROUTER_GLOBAL)`, the receivers over every node --
the number of shards by construction, one per receiver, and the same on every node, since the
Router refuses a team whose nodes differ in shape.  The hash buys the
producer a fan-out with no 64-bit division on the hot path; `murmur64` is a bijection, so distinct
images stay distinct.

**The shard.**  `DirectDict` is linear probing over 8-byte slots, `OCCUPIED | check << n | preimage`,
zero empty, `n <= 63`, a plain container that knows nothing of the problem, with four operations:
`insert`, `probe`, `prefetch` and `clear`.  The first three take the key **mixed once more** by the
dict thread, `h = murmur64(key)`: the run starts at `(h * n_slots) >> 64`, from `h`'s top bits, and
its tag -- bit 63 and `h`'s low `63 - n` bits, shifted above the preimage -- is what every slot
holding that key carries.  The second mix is required, not cosmetic: routing consumed the key's top
bits, so every key a shard receives lies in one interval of 2^64 / n_dicts, and a run cut from the
key itself would crowd one n_dicts-th of the shard.

`FILL`: `insert` into the first empty slot of the run; the run of a full shard never ends, which
`fill <= 0.9` rules out.  `PROBE`: `probe` fills a caller-owned `vector<u64>`, cleared first, with
the preimage of every slot of the run carrying the key's tag -- the **candidates**, each a collision
(`N_COLLISIONS`) bar a check-bit false positive, one slot in 2^(63-n).  The dict thread's `retire`
then tests each with the wrapper's `good(x, y)`, which is `is_good_pair(x, y) && f(x) == g(y)` --
for a collision problem `x != y && is_good_pair(x, y) && f(x) == f(y)`, `is_good_pair` being
symmetric by contract so the pair is judged in the order it came; and a collision search probes the
domain it has just filled, so every preimage of the chunk meets **itself**, `N_INSERT` self-matches
per phase that `x != y` rejects but `N_COLLISIONS` has already counted.  The predicate comes first
because it is cheap and rejects all but the golden pair; the evaluation, the guard of the answer
against the check bits' false positives, is spent only on what the predicate accepted.
`is_good_pair` is therefore **asked about pairs that do not collide**, which `problem.hpp` requires
it to answer.  What `good` accepts goes to `set_golden(x, y)` at once.  No collision queue and no
hand-off: a match costs about what a probe does, and a hand-off would cost more than it saves.  The
vector is the dict thread's own, one per round, sized once by its first long run and never by more
than the run.

A dict thread takes its points one at a time, with `Router_Pop`, which reads the block **where it
lies** and releases it on its last point: nothing is copied, and a block may hold fewer than
`block_points` points.  In `FILL` a point is inserted as it comes: the phase runs at the machine's
random-write pace and a prefetch buys it nothing, as measured on grdix.  In `PROBE` a point
popped is not retired at once: the first slot of its run is prefetched and the point waits in the
thread's **ring of `--prefetch` points** (`Options::prefetch`, default 8, at most 64, 0 probes each
point as it comes), retired when the point that many pops later arrives, so that a slot's DRAM miss
overlaps the pops and retirements behind it instead of being sat through alone -- on grdix, 192 dict
threads at one miss in flight each were the machine's pace.  The ring spans blocks and is
**drained**, every waiting point probed oldest first, once `Router_Pop` finds nothing and
`Router_Test_drained` is true, before the shard is cleared and before the thread returns: nothing is
dropped and the round's end is unchanged.  The ring is private to its dict thread, an array on its
stack.

### 1.4 The epilogue: one `MPI_Allgather` per phase

Thread 0, after `Router_Reset` and before the engine's barrier, builds its node's record and
exchanges it with every other node in a single `MPI_Allgather` of `REC_WORDS` `u64`:

| words | what |
|---|---|
| `0 .. N_COUNTERS-1` | the node's tallies, summed over its threads |
| `REC_ROUTER ..` | `ROUTER_STATS_SIZE` words: the node's `Router_Stats` for the phase |
| `REC_FOUND`, `REC_X`, `REC_Y` | 1 and the pair, if this node found a golden pair |

Every node then reads the same verdict out of the same records: **the lowest rank with `REC_FOUND`
provides the answer**, which thread 0 writes over the node's own golden pair: `found` and `golden` are
this node's until the epilogue and the verdict after it, the same on every node.  There is no stop
flag -- every thread leaves the loop when `found` is set or the domain is exhausted,
`phase == PROBE && round + 1 >= R`, a test it makes from its own loop.  Rank 0 also sums
the records and prints the round report, which is therefore exact.  Only thread 0 calls MPI, and it
does so while every other thread of the node waits at the engine's barrier; the Router sends nothing
on the communicator outside a round and `Router_Reset`'s barrier, so the `MPI_Allgather` cannot meet
one of its messages.

There is **no early exit inside a phase**: a solution found during a phase stops the search at the
end of that phase.  Cutting a phase short would need a control channel of its own, and it can save
at most one phase.

### 1.5 State, tallies and loss

| lives | where |
|---|---|
| a thread's handle, a dict thread's shard, a producer's buffers, thread 0's stats and record arrays | a local of that thread's scope, built by that thread |
| the tallies, the golden pair -- the node's own, then the verdict | `Shared`, one per rank, built before the team |

`Shared::tally[tid].ctr[]` is one `u64[N_COUNTERS]` per thread on a cache line of its own, **plain,
never atomic**: its owner writes it, thread 0 reads it after `Router_Reset`'s first team barrier,
which is what makes the writes visible.  The golden pair is the one exception, `set_golden` taking a
mutex and an `std::atomic` flag, because any dict thread of the node may find one at any moment.
Thread 0 then writes the verdict over it in the epilogue with no lock at all, every other thread of
the node being at the engine's barrier.

**Counters.**  `N_EVAL` (producers: evaluations, one point pushed each); `N_INSERT`, `N_PROBE`,
`N_COLLISIONS` (dict threads).  Everything about the communication -- points pushed and delivered,
bytes and messages on the wire, blocks stalled, the service thread's turns, the time the producers
spent waiting for a free block or closed, the dict threads without a block to read and the service
moving them --
is the Router's own tallies, in the same record; the engine keeps no counter of its own for it.

**Lossless.**  The Router loses nothing: a `Router_Push` may block, but no point is ever dropped.
A dropped point would be a missing entry or a missing probe, and "no solution" is a proof only when
none was dropped.  The round report prints the Router's two stall counters when they are nonzero,
and complains if the points pushed and the points delivered do not agree.

**The live line** is rank 0's own node, read from its tallies without synchronisation while the
phase runs, and scaled by the number of nodes: approximate on purpose.  The round report comes from
the `MPI_Allgather` and is exact.  Do not expect the two to agree to the unit.

**What limits the phase** is the Router's verdict, `Router_Bottleneck` over its stats (`router.3`):
the share of the phase the senders spent waiting for a free block or closed, the receivers with
nothing to pop and the service in turns that moved something, and one word -- `senders` when only
the receivers waited (the points are produced too slowly), `receivers` the other way round, `service/network`
when both waited (the blocks are in flight or in the service's hands; the busy share says which),
`balanced` when neither did -- by a 5 % threshold.  The live line ends with it, from rank 0's node;
the round report's `LIMITED BY` line has it from the summed records, exact.  The waits are timed at
their ends only, never per point, and the drain tail is receiver wait, so the verdict means something
for phases long against their closure; round 0's `FILL` says `receivers` while the dict threads still
touch their shards, which is true.

### 1.6 Answer

The golden slot holds `(x, y)`.  `claw_search` returns `(x, y)` with `f(x) == g(y)`,
`collision_search` a pair of distinct preimages of one image.  With no solution the search
returns nothing after `R` rounds, and that is a proof of absence -- unless `--nrounds` cut it short,
which the last line says.  The demos' exit status is that outcome.

## 2. PCS, on the Router

`mitm::pcs` (`include/pcs/`): van Oorschot-Wiener parallel collision search over a
distributed dictionary of trail endpoints.  The problem becomes one random function of its
range onto itself, `x -> f(mix(i, x))`; a **round is one version `i` of that function**, and
it ends when `beta * w` distinguished points have been found over **all** the nodes
(`--beta`, default 8).  The next round draws a fresh `i`, so a run is a sequence of
independent functions and the pair is found in one of them.  Unlike §1 this is a Monte-Carlo
search: it stops when it finds the pair, and "not found" is never a proof.

### 2.1 The team

Exactly §1.1's shape: one OpenMP team of `1 + I + W` threads per rank, created once for the
whole run, every thread calling
`Router_Init(role, ROUTER_GROUP_AUTO, comm, ROUTER_TAG, &opts.router)` once and keeping its
handle.

| thread | Router role | does |
|---|---|---|
| 0 | service | all the MPI: `Router_Progress` until `Router_Test_quiescent`, the control channel, the epilogue, the printing |
| `1..I` | receiver | the dict thread: owns one shard, probes what the Router delivers and hands the hits to its walkers |
| `I+1..I+W` | sender | the walker: walks `vlen` trails, pushes every distinguished point, and resolves the hits its dict thread found |

**The Router owns placement**, exactly as in §1: PCS pins nothing and never touches hwloc.  It
only *reads* the plan, through `Router_group` and `Router_num_groups`, and **the team's shape**
through `Router_size(rt, role, scope)` and `Router_rank(rt, scope)`: how many walkers, dict
threads or nodes there are (`ROUTER_NODE` or `ROUTER_GLOBAL`), and a thread's own rank among
those of its role.  PCS stores none of it -- `pcs::Params` is the dictionary's geometry and the
round's quota, nothing about the team -- and asks wherever it needs one: the chain-index stride
and a walker's first chain (§2.3), the routing modulus (§2.3), the pass below, the size of the
epilogue's record array (§2.5) and the controller's end-of-round to every node (§2.2).

**Every per-thread object is built by its owner right after `Router_Init`**, i.e. once pinned,
so the zero-fill is the first touch of every page: a dict thread's shard and its collision
queue, a walker's HyperLogLog registers.  A walker's lanes are
the round's own, built afresh at its start.  That
costs **one team barrier** §1 does not need, after those allocations and before the first
round, because the groups and the queues have to be published to the team.

**Which dict thread a walker resolves for.**  A hit is handed over, not resolved on the spot
(§2.4), so every walker is paired with exactly one dict thread and every dict thread needs at
least one walker.  The Router's groups do not give that for free: a group is senders and
receivers in one cache domain, so it may hold several receivers, or none, and the Router exposes
no query for "the receivers of a group" -- only `Router_group`, one thread's own, and
`Router_num_groups`.  So each worker publishes its group into the `Shared`, at its local rank
(`Router_rank(rt, ROUTER_NODE)`; the dict threads' slots first, then the walkers'), and after the
startup barrier every thread runs the same deterministic pass over those groups:

1. **the groups that hold a receiver**, in group order: each of that group's senders, in
   order, goes to the **poorest receiver of that group** -- fewest walkers so far, ties to the
   lowest index.  Serving a group from inside itself is the point of the pass: where a group
   is a cache domain, the queue and the two threads touching it stay in one cache;
2. **the senders left over**, whose group holds no receiver at all, in order, each to the
   **poorest receiver of the node**;
3. and if a receiver is *still* without a walker, the run **stops**.  There is no third pass
   and no borrowing: such a shard would find collisions nobody resolves and drop every one of
   them, so thread 0 reports which knobs to turn and calls `MPI_Abort`.  `I > W` cannot work
   at all and is refused before the team is even built.

### 2.2 The round, and the control channel

A round ends on a count no node holds, so PCS keeps what §1 does without: a **controller on
rank 0**, driven by that rank's service thread between two Router turns, and three tags of its
own on the Router's communicator.  Every message is `MPI_Bsend` into a buffer the `Control`
attaches when thread 0 builds it, right after `Router_Init`, and gives back when the run ends,
so a report never waits on the controller; all of it is thread 0's -- the same thread the Router
does its MPI on, which is what keeps `MPI_THREAD_FUNNELED` legal.

| tag | direction | payload | when |
|---|---|---|---|
| `TAG_REPORT` | node -> rank 0 | `REC_FOUND` words: the node's tallies then its `Router_Stats`, as a **delta** since its last report | when `N_DP` has advanced by `report_points`, or `ping_delay` has elapsed -- tested every 256th service turn |
| `TAG_END_ROUND` | rank 0 -> every node, itself included | zero length | once per round, when `reported[N_DP] >= beta * w`, or a node has said it holds the pair |
| `TAG_SOLUTION` | node -> rank 0 | the version and the two points | once, when a walker sets the node's golden pair |

`report_points` is `beta * w / (64 * nodes)`, so a node owes about 64 reports a round beside
the ones the clock asks for: on a timer alone a short round overshoots its quota, and on
volume alone a slow one goes quiet.  Rank 0 keeps one `MPI_Irecv(MPI_ANY_SOURCE)` posted for
each of the two inbound tags and **drains both to empty** on every control turn, which is what
keeps a report matched late in its own drain from being credited to the next round; the
`reported[]` it accumulates is approximate on purpose and decides one thing only, when the
round closes.  Every node, rank 0 included, keeps one `MPI_Irecv` posted for
`TAG_END_ROUND` and reposts it the moment it matches, which is always before `Router_Reset`,
so the next round's token can never overtake it.  `TAG_SOLUTION` only ends the round early --
the pair itself reaches every node through the epilogue.

One round, on each thread:

1. zero its own tallies;
2. its role's work:
   - a **walker** runs its `vlen` lanes: where an instruction ends the lane takes the next one,
     a candidate its dict thread queued if there is one and a fresh trail otherwise, and every
     distinguished point is pushed as it is found.  It tests `round_over` every cycle; when it
     is up it calls `Router_Close` **once** and then runs the machine dry;
   - a **dict thread** takes points one at a time with `Router_Pop`, probes each into its
     shard and fills a run of candidates, handing the run over as soon as there is nothing
     left to probe.  It loops until `Router_Pop` finds nothing and `Router_Test_drained` is
     true, then marks its queue final and empties its shard;
   - **thread 0** turns the Router until the node is quiescent, serving the control channel
     every 256th turn -- which is how the end of the round arrives and the walkers get to
     close -- and then reads `Router_Stats`, which `Router_Reset` would clear;
3. `Router_Reset`, which every thread of every node calls;
4. thread 0 alone: the epilogue (§2.5);
5. `#pragma omp barrier`, the engine's own and only one: it publishes the verdict and the next
   round's header, and lets the next round zero the tallies thread 0 has just read;
6. every thread leaves the loop if `found` is set or `--nrounds` versions have run -- a test it
   makes for itself -- or starts the next round.

**The drain, and why it needs a handshake.**  A candidate walked in the next round would be
walked under the wrong version of the function, so **every candidate must be retired before
`Router_Reset`**.  A dict thread keeps producing candidates until `Router_Test_drained`, which
is node-wide and needs every sender of every node closed; so each dict thread carries one
`done` flag, set on release when its pop loop ends, and a walker's drain ends when **its dict
thread is done, the queue it shares with its fellow walkers is empty and its own hand is
empty** -- `done` read before the queue, or a run pushed between the two tests would be left
behind.  That flag and `round_over` are the whole of it: the previous core's
`RUNNING / HOLD / HELD / DRAIN / QUIESCENT` state machine is gone.  `done` is set *before* the
shard is emptied, so a walker's drain overlaps that flush.

**`Router_Close` comes first, not last**, and that ordering is forced: the round cannot end
until every sender has closed, so `Router_Test_drained` -- and with it `done` -- is a flag the
walker's own silence has to earn.  A walker that waited for its queue before closing would wait
for ever, and would keep feeding the dictionary the very candidates it was waiting to stop
seeing.  Closed, it still runs every lane: the trails keep the machine full and their points go
nowhere, being neither tallied nor pushed.

**Then the flush.**  Once the queue is final and empty, nothing new can reach a lane, and the
walker steps `4 * dp_max_it` more cycles -- inside which every collision still on a lane must
have run its course, since each unknown length is re-walked in at most `dp_max_it` and the
align and the march together cost `len0 + len1`, the march two cycles a step.  It then asserts
that every lane is back on a trail.  A bounded count and an assertion, rather than a live count
of busy lanes: the machine's one invariant is that a lane always has an instruction, and asking
it how many are collisions would put back the scan §2.4 exists to avoid.

**The round header travels on no wire.**  Every rank holds the same `PRNG` -- the driver draws
one seed and broadcasts it -- so thread 0 draws `Header{i, root_seed}` for the first round
before the team and every later one in the epilogue, and the engine's barrier publishes it.
The verdict is identical on every node too (§2.5), so the PRNGs stay in step for ever and
`i` needs no broadcast.

### 2.3 Points and the shard

A point is the two words the Router carries: `key` is the trail's endpoint in full, `val` its
chain index in the low `jbits` with the trail's length above that, saturating -- the largest
value means "at least this long, exact length unknown".  A chain index is all a start point
needs, since chain `j` starts at `root_seed + j * multiplier`: nothing about a trail's start
travels, and any walker on any node can reconstruct it from the header.

**Routing.**  `dest = key % n_dicts`, where `n_dicts` is `nodes * I` -- not a Router query, the
same shard count every node derives from its own `Options`.  The shard then keys on
`key / n_dicts`, which is exactly the information the modulo left it, so no hash is needed
anywhere: a trail endpoint is already a function value.

**The shard.**  `PcsDict` is one slot per endpoint, direct-mapped, no probing: `key / n_slots`
above the length, the length in 8 bits above the chain index.  A hit is only ever a hint --
the walker re-walks both trails anyway -- so a slot is simply overwritten, and only by a trail
**at least as long** as the one it holds.  A length that saturated on the wire is stored as
unknown whatever the 8-bit field could have held, so the two saturations stay consistent and
the dictionary never reports an exact length it does not have.  The shard is emptied at the
end of every round.

### 2.4 The collision queues

`CollisionCandidate` is what a dict thread hands to its walkers on a hit:

| field | meaning |
|---|---|
| `i` | function version; the walker asserts it is the round's |
| `seed0`, `len0_maybe` | the incoming point's chain index and trail length; `len0_maybe == 0` means the length had saturated **on the wire** |
| `end` | the shared endpoint, in full: what a re-walked trail must reach |
| `seed1`, `len1_maybe` | the point already in the slot; `len1_maybe == 0` means the stored length had saturated (8 bits) |

**Either length can be unknown, or both**, and `0` is free as the marker because a trail's
length is at least 1.  The walker recovers an unknown length by re-walking that trail to its
distinguished point, which must be `end` and must be reached within `dp_max_it` steps;
`N_MEASURE` counts the re-walks, `BAD_DP` and `BAD_WALK_NONCOLLIDING` the two ways one fails.

**One queue per dict thread, and candidates cross it in runs.**  Collisions are not rare: at
`beta = 8` a fair share of every probe hits, so the queue carries a fixed share of the traffic
the shards do.  One queue for the whole node, pushed by every dict thread, is what that share
cannot survive -- so there is one per dict thread, its consumers are the walkers §2.1 gave it,
and both sides move up to 64 candidates per lock acquisition.  The mutex stays: what made a
single queue costly was that there was one of it, and that its `count` line was written
millions of times a second, so every walker's poll missed.  Batched and split, that line is
written thousands of times a second and a poll finds it valid; the emptiness probe a walker
runs before taking the lock is a relaxed atomic read, so an idle queue costs no lock traffic
at all.

**The queue is the one lossy thing in the engine, and that is forced.**  A full queue
truncates the run and the dict thread tallies the rest as `DROP_COLL`.  It cannot block
instead: `Router_Push` spins waiting for a free block, so a walker inside a push is draining
nothing, and a dict thread blocked on a queue whose only consumer was that walker would
deadlock the pair.  A dropped candidate costs one collision, which a Monte-Carlo search can
afford; the round report prints the count whenever it is not zero.  Nothing else is lost: the
round's end retires every candidate that reached a queue (§2.2).

**How a walker retires them.**  A round spends about `2/beta` of its evaluations locating
collisions, and resolving one point at a time makes each of those cost `vlen` times what a
walked one does -- most of a walker's time, once the dictionary fills.  So a walker is **one
machine of `vlen` lanes** that walks trails and resolves collisions on the same lanes: one
`vmixf` per cycle advances every lane, a lane runs one instruction to its end -- a fresh trail
to its distinguished point, or the resolution of one candidate -- and a lane that ends one is
redispatched at once, **a pending collision before a fresh trail**.  Nothing is preempted, and
there is no batch: trails and collisions end at different times, so the lanes desynchronise by
themselves.  A collision holds **one lane for its whole life** and walks its two chains through
three phases on it:

| phase | what it does |
|---|---|
| `MEASURE0`, `MEASURE1` | a length is unknown: re-walk that chain to its distinguished point to learn it -- chain 0, then chain 1 if both are.  A chain gives up after `dp_max_it` steps (`BAD_DP`), and abandons the candidate if it lands on a point other than `end` (`BAD_WALK_NONCOLLIDING`) |
| `ALIGN` | rewind both chains and step the longer one until both are the same distance from their shared endpoint; the shorter waits at its start, in the lane's registers |
| `MARCH` | step both and compare, until they meet (the collision) or the shorter trail runs out (`BAD_WALK_NONCOLLIDING`), **alternating the two chains on the lane**: chain 0 on the even cycles, its image saved, chain 1 on the odd ones, then the comparison.  Equal positions on entry mean one trail is a suffix of the other (`BAD_WALK_ROBINHOOD`) |

Alternating costs exactly what two lanes would -- two lane-cycles per step -- and only doubles
the latency of that phase; in exchange no instruction ever needs a second lane, so there is
nothing to allocate, nothing to wait for and no way to deadlock.  Every phase advances one chain
by one evaluation per cycle, and a lane that ends an instruction takes the next one on the spot,
so **no lane is ever idle**: a cycle is always `vlen` evaluations, `N_EVAL` counts them as such,
and the machine keeps neither a free list nor a live count -- which is what the round's end is
shaped around (§2.2).  A problem with `vlen == 1` is the same machine with one lane.

**Dispatch is collisions first, absolutely**: while candidates are pending -- the walker's
private run, `max(64, vlen)` deep and topped up from the queue once a cycle whenever fewer than
`vlen` are in hand, so that a cycle ending every lane's instruction still finds one for each --
every freed lane takes one, so a walker stops producing while its queue is deep, and starts a
trail on every lane its collisions do not need.  A saturated length is re-walked, not recorded: a record would cost
`dp_max_it + 1` words per lane, which bitsliced DES's 512 lanes cannot afford, to save one
re-walk per saturated length -- about none per collision up to `1/theta ~ 26`, 0.85 at
`1/theta ~ 212`.

### 2.5 The epilogue: one `MPI_Allgather` per round

Thread 0, after `Router_Reset` and before the engine's barrier, builds its node's record and
exchanges it with every other node in a single `MPI_Allgather` of `REC_WORDS` `u64`.  Its
first `REC_FOUND` words are the same layout a progress report carries, which is why one
function builds both:

| words | what |
|---|---|
| `0 .. N_COUNTERS-1` | the node's tallies, summed over its threads |
| `REC_ROUTER ..` | `ROUTER_STATS_SIZE` words: the node's `Router_Stats` for the round |
| `REC_FOUND`, `REC_I`, `REC_X0`, `REC_X1` | 1, the version and the pair, if this node found one |

Every node then reads the same verdict out of the same records: **the lowest rank with
`REC_FOUND` set provides the answer**, which thread 0 writes over the node's own golden pair:
`found` and `golden` are this node's until the epilogue and the verdict after it, the same on
every node.  The flip never becomes a `TAG_SOLUTION`, because no control turn follows an
epilogue that sets `found`.  There is no stop flag -- every thread leaves the loop when `found`
is set or `nround >= max_versions` (`--nrounds`), a test it makes for itself.  Rank 0 also sums
the records and prints the round report, which is therefore
exact -- the live line, read from the reports, is not.  Then the walkers' HyperLogLog
registers, merged and zeroed by thread 0 and `MPI_Reduce(MPI_MAX)`'d to rank 0: 64 KB a round
buys the round's count of **distinct** collisions, which is the number that says whether the
parameters are any good, and how a search that finds collisions but never the pair is told
from one that finds the same collision over and over.

Only thread 0 calls MPI, and it does so while every other thread of the node waits at the
engine's barrier.  The Router sends nothing on the communicator outside a round and
`Router_Reset`'s barrier, so neither collective can meet one of its messages; a peer's control
message can arrive during them, and waits in MPI's own queues until rank 0's next drain.

There is **no early exit inside a round**: a pair found while the round runs ends it through
`TAG_SOLUTION`, at the next control turn, not at once.

### 2.6 State, tallies and loss

| lives | where |
|---|---|
| a thread's handle, a dict thread's shard and queue, a walker's registers and lanes, thread 0's stats, records and controller | a local of that thread's scope, built by that thread |
| the tallies, the queues' pointers, the published groups, the header, the golden pair -- the node's own, then the verdict | the one `Shared`, built before the team |

`Shared::tally[tid].ctr[]` is one `u64[N_COUNTERS]` per thread on a cache line of its own,
**plain, never atomic**: its owner writes it, thread 0 reads it after `Router_Reset`'s first
team barrier, which is what makes the writes visible -- and reads it *without*
synchronisation during the round, for the report and the live line, which is why both are
approximate.  The exceptions are the golden pair (a mutex and an atomic flag, because any
walker may find one at any moment; thread 0 then writes the verdict over it in the epilogue
with no lock at all, every other thread being at the engine's barrier) and the two flags on
the hot path, `round_over` and each dict thread's `done`.

**Counters.**  Walkers: `N_EVAL`, `N_DP`, `N_POINTS_TRAILS`, `N_COLLISIONS`,
`COLLIDING_LEN_MIN`, `COLLIDING_LEN_MAX`, `N_MEASURE`, `BAD_DP`, `BAD_COLLISION`,
`BAD_WALK_ROBINHOOD`, `BAD_WALK_NONCOLLIDING`.  Dict threads: `N_PROBE`, `BAD_PROBE`,
`DROP_COLL`.  Everything about the communication is the Router's own tallies, in the same
record -- the walkers' and the dict threads' waits and the service's busy time included -- and
the engine keeps no counter of its own for it: the previous core's three drop counters have no
meaning on a transport that drops nothing.

**Loss.**  The Router loses nothing, and the round report says so -- it prints the Router's two
stall counters when they are nonzero, and complains if the points pushed and the points
delivered disagree.  A round closes on points **found**, not on points inserted, so a run
whose collision queues overflow burns rounds against a nearly empty dictionary and still
reports them complete: the report's `ROUTED` line -- how many of the points found reached a
shard -- is the number to read.

**What limits the round** is the same verdict as the direct engine's (§1.5), `Router_Bottleneck`
over the Router's stats: on the live line from the reports, so over every node but lagging them;
on the round report's `LIMITED BY` line from the summed records, exact.  A walker's wait is its
`Router_Push` and its time closed -- the time it spends resolving collisions is its own work, not a
wait -- so a round whose walkers resolve more than they walk reads as `senders`, which is right.

### 2.7 Answer

The golden slot holds `(i, x0, x1)`: the version, and the two inputs of one evaluation of the
mixed function.  `unmix` turns them back into the problem's own terms, putting `f`'s side
first for a claw, and the entry points re-verify the answer against the raw problem outside
the engine.  With no solution the search returns nothing after `--nrounds` versions, and that
is **not** a proof of absence -- it is where the direct engine's §1.6 differs.  The demos'
exit status is that outcome.
