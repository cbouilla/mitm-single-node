# Routing is the bottleneck on a large node

Measured on `grdix-4` (2 sockets, 256 cores / 512 PUs, 2 NUMA nodes of 256 PUs each,
1 TB) in September 2026, on `double_speck64_demo`.  The L3 layout below is inferred
from the placement the banner reports, not from `lstopo`; check it before relying on
it.  Everything below is from
`--n 40 --ram 4G`, so `w = 5.0e8` slots (`2^28.90`) and the auto-tuned trail length is
`1/theta == 18.76`.

**The engine delivers about 4 M distinguished points per second to its dictionary, on
a machine whose walkers can produce 760 M/s.  Nothing tried so far moved that number,
because every point on a node passes through one thread.**

## 1. What was measured

`-O3 -DNDEBUG` throughout; the first `-O0` run is left out, its rates mean nothing.
"offered" is `walkers * (f/s) / L`; "absorbed" is `inserters * (probe/s)`, both as the
round report prints them.

| walkers / inserters | threads | offered | walker-queue drops | absorbed | dict occupancy | collisions | coll/s |
|---|---|---|---|---|---|---|---|
| 507 / 4 | 512 | 757 M/s | 95.9% | 3.9 M/s | 3.1% | 2^16.19 | 13.8 K |
| 250 / 4 | 255 | 456 M/s | 92.4% | 4.4 M/s | 7.3% | 2^17.88 | 27.0 K |
| 24 / 32 | 57 | 48 M/s | 90.4% | 4.6 M/s | 76.7% | 2^24.17 | 219 K |
| 4 / 4 | 9 | 2.8 M/s | 0.88% | 2.8 M/s | 8x overfilled | 2^29.33 | **477 K** |

Two readings that are easy to get wrong:

- **Wall-clock per round is not the metric.**  A round closes on distinguished points
  *found* (`controller.hpp:96`, against `N_DP`), and in the first three rows over 90%
  of those were dropped before reaching a shard.  The 512-thread run finished a round
  in 5.4 s having inserted 3% of `beta*w`; the 9-thread run took 1436 s and filled the
  dictionary eight times over.  Collisions per second is the objective, and on it the
  9-thread configuration beats the 512-thread one by **34x**.
- **The comm thread's apparent 31-34 M DP/s in rows 1 and 2 was 87% *failed* pushes** --
  a `head.load` and a return, no write -- because the inserter queues were permanently
  full.  Its real delivery has been 2.8-4.6 M DP/s in every configuration ever run.

## 2. Why this is the whole game

```
rounds        = 1.8 * 2^n / w
DPs per round = beta * w
total DPs     = 1.8 * beta * 2^n          <- independent of w, theta, alpha
runtime       = 1.8 * beta * 2^n / D_total
```

The number of distinguished points the whole job must route is fixed by `n` alone.
**Routing throughput is the attack's speed, one for one.**  At `--n 40` that is
`1.58e13` points: 41 days at 4.5 M/s, 5.8 hours at the 760 M/s the walkers already
produce.

Three consequences:

- **Difficulty cannot buy anything.**  `total f-evals = 1.8 * beta * 2^n * L` with
  `L = sqrt(2^n / w) / alpha`, so lowering `theta` cuts the DP rate and raises the
  f-work by the same factor.  `runtime * DP rate` is a constant.
- **Memory makes it worse.**  DP work is flat in `w`; f-work falls as `1/sqrt(w)`.  A
  node with 1 TB and 256 cores is by construction one whose walkers are nearly free
  and whose router is everything.  Using the full 1 TB instead of 4 GB multiplies the
  DP rate by `sqrt(250) = 16`.
- **Scale makes it worse too.**  `w = ram_per_node * N / 8`, so
  `D_node = R_f * alpha * sqrt(w / 2^n)` grows as `sqrt(N)`.  The funnel tightens as
  the machine gets bigger, in both directions at once.

Where one comm thread at `C = 4.5 M DP/s` is enough, with `alpha = 2.5` and
`R_f = 1.4e10 f/s` per node (507 walkers x 28 M, SMT on, dictionary still empty), from
`w <= 2^n * (C / (alpha * R_f))^2`:

| RAM/node | dictionary usable from | nodes before the router saturates, at n = 64 |
|---|---|---|
| 4 GB | n >= 55 | ~600 |
| 1 TB | n >= 63 | **~2** |

That last cell is the finding.  With the memory this machine actually has, the
one-comm-thread design runs out at about two nodes.

## 3. Where the comm thread's time goes

`route_walker_queues` (`engine.hpp:82`) sweeps every walker queue with `pop_bulk(64)`,
then `deliver_local` (`engine.hpp:74`) writes each 24-byte point into the ring of the
shard that owns it.  Measured cost per *delivered* point, 32 inserters: **220 ns**.

- **Distance.**  One thread pulls from up to 507 queues spread over ~16 L3 domains and
  both sockets, then scatters into rings whose consumers are equally spread.  Nearly
  every line it touches arrives from another die.
- **True sharing on the ring data.**  With 32 inserters the queues are never full
  (`DROP_INSERTERQ == 0`), so the consumer sits right behind each write; 24-byte points
  put 2-3 of them per 64-byte line, and half the consumers are across the socket.
- **One `tail.store(release)` per point** (`spsc.hpp:44`), which each idle inserter
  re-reads from its `in.empty()` spin.  `SPSCQueue` has `pop_bulk` and no `push_bulk`;
  that asymmetry is the bug.
- **Two runtime 64-bit divisions per point** in `deliver_local` -- `x % n_nodes` and
  `x / n_nodes % inserters_per_node` -- with `n_nodes == 1`, because both are runtime
  members.  The inserter adds three more (`inserter.hpp:53`, `:54`, `:101`).

All but the last scale *with* the inserter count, which is why 4 -> 32 inserters made
the comm thread 7.5x slower instead of faster.

Two smaller effects, for completeness.  Thread placement (`parameters.hpp:143`) is
NUMA-aware and SMT-blind: it takes the next free CPU of a NUMA node, so once it passes
that node's primary threads it hands walkers the siblings of the comm thread and of
every inserter.  Measured by leaving the siblings idle (`--walkers-per-node 250`):
worth ~10% on those threads, while SMT is worth +66% to the walkers -- so the fix is
SMT-aware placement, not fewer threads.  And an inserter probe costs ~900 ns against a 1 GB shard, against ~150 ns
for one dependent random access with a page walk; huge pages and software pipelining
(the TODO at `inserter.hpp:82`) are untested on this hardware.

## 4. What is not the answer

- **More inserters.**  Two targets, both unreachable.  Matching what the walkers
  produce needs `I/W = alpha * (R_f/R_p) * sqrt(w / 2^n) = 1.66` -- more inserters than
  walkers, 158 of the 254 workers.  Matching only what the comm thread delivers needs
  `31M / 1.1M = 28`.  Going from 4 to 32 made things worse, per section 3: the balance
  is right, the queue implementation punishes it.
- **Fewer walkers.**  The 9-thread run wins only because it stops overwhelming the
  router.  It is a diagnosis, not a configuration to ship.
- **Walkers delivering into the shards themselves** (MPSC, bypassing the comm thread).
  That removes the `1/N` local fraction and leaves `(N-1)/N` exactly where it was.  A
  single-node trick, worthless at scale, not worth a protocol change.
- **Turning `alpha` down.**  Section 2: it trades one for one.

## 5. Ideas worth trying

### A. One MPI rank per L3 domain -- zero code

`R` ranks per node is `R` comm threads, and the routing work splits `R` ways; whether
a point's shard is on this node stops mattering.  Binding a rank to an L3 domain makes
its walkers, its staging buffers and its shards all domain-local, so the per-router
rate should rise as well as the count.

```bash
mpirun -np 16 --map-by l3cache --bind-to l3cache \
    build/examples/double_speck64_demo --n 40 --ram 250M \
    --inserters-per-node 4 --walkers-per-node 8
```

`--ram 250M x 16` keeps `w` at the same 4 GB total, so the insertion rate is directly
comparable.  Check the banner reports 32 CPUs in the affinity mask, not 512: nothing
coordinates placement across ranks, so ranks that inherit the same mask will silently
pin to the same CPUs.

Cost: send buffers grow as `O(ranks^2)` (16 ranks x 1000 nodes is 1.15 GB per rank),
a round barrier over `16N` ranks instead of `N`, and `16N` progress reports converging
on rank 0.  Fine to ~10-20k ranks.  This contradicts the "one MPI rank per node"
convention in `CLAUDE.md`, which is exactly what makes it interesting.

### B. Router threads, one per L3 domain -- the in-process version of A

A router owns the walker queues, the staging buffers and the inserters of its domain;
thread 0 keeps every MPI call (`THREAD_FUNNELED` is why the funnel exists) but does
**no per-point work**: it posts an `Isend` for a full buffer, and hands an arriving
buffer to the router that owns it.

Two things this has to get right:

- **The destination key must name a single consumer.**  Today `dst = x % n_nodes`
  (`engine.hpp:92`), so an arriving buffer holds points for every shard on the target
  node and someone must scan and scatter it -- `O(D)` on one thread, on the *receiving*
  side too, which is half the problem and would survive the change untouched.  Key on
  `(node, domain)` or on the global shard index instead.
- **Buffer memory sets `buffer_capacity`.**  `routers * destinations * capacity`, and
  it wants to be L3-resident: ~32 MB per router, so ~850 points per buffer at 1600
  global shards.  At 16k shards that falls to 85, too small to amortise a message, and
  the key has to become two-level (`node`, then domain within it) to keep the
  destination count at `N * R`.

Over A this buys: no shared-memory copy for intra-node traffic, one dictionary and one
round barrier per node instead of 16, and `O(R^2 N)` buffers instead of `O((16N)^2)`.
It costs a real `PROTOCOL.md` rewrite -- a `ROUTER` role in §1, new channels in §3,
the `TAG_POINTS` addressing and wire layout in §2, and a new stage in the §4.5 drain
(routers must quiesce after the walkers and before the buffers flush).

**Measure A first.**  It is the same architecture implemented by the MPI runtime, and
it bounds what B is worth before anyone spends a protocol rewrite on it.

### C. Constants, worth having either way

Each multiplies whatever parallel routing gives.

- `push_bulk` on `SPSCQueue`: bucket a whole sweep by destination and push each bucket
  once -- ~48 points per `tail.store` instead of one, and contiguous writes per ring.
- Multiply-shift reciprocals for the five runtime divisions of section 3; all five
  divisors are fixed at construction.
- `MADV_HUGEPAGE` on the shard (which means replacing `vector<u64> A`, since its
  zero-fill is also the NUMA first touch), and the software-pipelined probe.

## 6. Accounting that hides all this

- **A round closes on points found, not inserted** (`controller.hpp:96`).  A starved
  run burns rounds against an almost-empty dictionary and reports them complete.  At
  the least the controller should say so when `DROP_WALKERQ` is a large fraction of
  `N_DP`.
- **`controller.hpp:179` divides `BAD_PROBE` by `N_DP`** -- points *found* -- when it
  is an inserter-side counter whose denominator is `N_PROBE`.  With 96% of points
  dropped the two differ by ~250x, and a printed "0.52% probe failure" was really
  ~100%.  The other four percentages on that line are walker-side, so `ndp` is right
  for them.
- There is still **no way to measure a shard probe in isolation**: `benchmark.hpp`
  covers `f`/`g`/`vfg` only, so `R_p` has only ever been seen through the full
  pipeline, tangled up with the queues and the router.

## 7. Numbers to re-measure

The resolver is now batched (`PROTOCOL.md` §3.2), so walkers stay near their
empty-dictionary rate much further into a round than they did when these runs were
taken -- which means they will offer *more* points and reach the router's ceiling
sooner, not later.  Everything in section 1 should be re-run before it is quoted.

Worth having per node, once: `double_speck64_bench` (the scalar/vector ratio -- 5.9x
on a laptop), the comm thread's delivered rate at 1, 4 and 16 inserters, and a probe
microbenchmark once one exists.
