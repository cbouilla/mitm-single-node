# Router benchmark: where the time goes

Two profiling sessions, both 2026-09-07.  **Session 1** (below) is a 32-core Intel
**grvingt** node; **session 2** (at the end) is a 256-core AMD **grdix** node, which
reproduces the analysis at 8x the core count and *corrects* several conclusions of
session 1 -- read both before trusting a number.  Where they disagree, the disagreement
is called out explicitly in session 2 §6.

## Session 1 -- grvingt (32 cores, Intel)

Profiling session, 2026-09-07, on a **grvingt** node (Grid'5000, Nancy): 2 sockets x
Intel Xeon Gold 6130 (16 cores / 32 PU per socket, 64 PU total), 2 NUMA nodes,
44 MB L3 (2 x 22 MB), Open MPI 5.0.7.  Binaries as built by
`cmake -DCMAKE_BUILD_TYPE=release` (`-O3 -DNDEBUG -march=native -ggdb`).  `perf` needed
`perf_event_paranoid = -1` (`sudo-g5k`).

Everything below is `router_bench`.  Some numbers come from out-of-tree variants of
`examples/router_bench.cpp` (noted where used); the library was not modified.

## Summary

Two very different regimes, and they must not be confused:

| configuration | sender loop | routed | ns/point | senders in `cpu_relax` |
|---|---|---|---|---|
| np 1, 25 senders + 6 receivers, both sockets | stock | 1.2 G pts/s | 20.9 | 2 % |
| np 1, 25 + 6, both sockets | no divider | 1.4 G pts/s | 18.4 | - |
| np 2 (one rank per NUMA node), 12 + 3, node-local destinations | no divider | **2.0 G pts/s** | **11.9** | ~0 % |
| np 2 (one rank per NUMA node), 12 + 3, uniform destinations | no divider | 0.68 G pts/s | 35.2 | **43 %** |
| np 2, 12 + 3, uniform destinations | stock | 0.69 G pts/s | 34.7 | - |

"stock" is `examples/router_bench.cpp` as committed; "no divider" replaces its per-point
`x % F` with a multiply-shift (see §3, which is why the two must be kept apart).

The local path is fast.  The moment half the points cross MPI, throughput drops 3x and
the senders spend nearly half their cycles spinning.  One rank per NUMA node with
node-local traffic is the best configuration measured, and it beats spreading one rank
over both sockets (2.0 G vs 1.4 G with the same sender loop).

## 1. The senders' `cpu_relax` is the headline (np 2, uniform destinations)

`perf record -F 999`, rank 0, steady state, samples attributed per instruction and the
three `pause` sites resolved with `addr2line -i`.  Of the 12 sender threads' cycles:

| where | share of sender cycles |
|---|---|
| **`router.hpp:837` — `stage_line`'s "a sealer is installing" wait** | **28.0 %** |
| **`router.hpp:769` — `take_fresh()`, the free-block stack is empty** | **13.5 %** |
| the atomic load in that wait loop | 1.7 % |
| `router.hpp:845-852`, the body of `Router_Push` | ~14 % |
| `x % F` in the benchmark (see §3) | 12.8 % |
| PRNG | ~13 % |

The mechanism, in order:

1. A push fills the L-th line of a destination's block, so that sender becomes the
   **sealer**: `stage_line` calls `take_fresh()`.
2. The free-block stack is empty, because blocks destined for the other rank are not
   returned until their MPI send completes.  The sealer spins at `router.hpp:769`.
3. While it spins, the destination's `next` word is left in the "installing" state
   (`k > L`), so **every other sender pushing to that destination convoys** behind it at
   `router.hpp:837`.

That is why 13.5 % of spin in `take_fresh` (one sender at a time) amplifies into 28 % at
the install wait across the other eleven.  It is head-of-line blocking through a single
global free stack: back-pressure from the remote path stalls senders whose points were
going to a *local* receiver.

For contrast, the same measurement at np 1 (`n_nodes == 1`, no MPI at all) puts senders at
**2.1 %** `cpu_relax`.  The spin is entirely a consequence of the inter-rank path.

## 2. The inter-rank path: one saturated service thread doing cold cross-socket copies

Measured 2.70-2.84 GB/s of inter-rank traffic, 41-43 K messages/s of 65600 bytes.

The reported "service 13-32 % busy" is **misleading** — it is a fraction of *turns*, not of
cycles.  A user-space-only profile also makes the service thread look idle (906 M user
cycles in 4.8 s).  With kernel symbols included it is **saturated**: 9.7 G cycles in 4 s
(~100 % of a core), of which

- **73.5 % `rep_movs_alternative`** (the kernel's memcpy), reached through
  `process_vm_rw_core` — Open MPI's `sm` BTL using CMA single-copy;
- the rest largely per-transfer overhead: `__get_user_pages`, `unpin_user_pages`,
  `__ptrace_may_access`, `apparmor_ptrace_access_check`.

Corroborated by wall-clock: `sys` time is **5.28 s** with MPI traffic and **0.15 s** with
node-local destinations only.

**This is not the router's MPI usage.**  A standalone MPI streaming test — one thread per
rank, bidirectional, 64 KB messages, a send window of 4 and 32 posted receives, i.e. the
router's own pattern, rotating through N distinct send buffers — gives:

| distinct send buffers rotated | GB/s per rank |
|---|---|
| 4 (hot in cache) | 6.09 |
| 32 | 4.30 |
| 256 | 3.83 |
| **725 (= the router's pool)** | **2.74** |

The router's 2.7 GB/s *is* the achievable rate for a single thread copying cold,
cross-socket 64 KB blocks.  Two structural causes, both outside the MPI calls:

- all inter-rank copying is funnelled through one thread (`MPI_THREAD_FUNNELED`);
- the block pool (725 blocks x 65600 B = 47.6 MB) is far larger than a socket's 22 MB L3,
  so every copy reads and writes cold DRAM, and across UPI.

Setting `--mca btl_sm_single_copy_mechanism none` changes nothing (34.9 vs 35.3 ns/point),
which is consistent: the limit is the copy itself, not the mechanism.

### Knobs that matter

- **`--credit 4` is too tight.**  16 or 64 gives +13 % (34.7 -> 30.7 ns/point).  Beyond 16
  there is nothing more.
- **Bigger blocks help, small blocks are ruinous** (np 2, uniform, `--swc` pinned at 512 so
  the two are separated):

  | `--block` | routed | ns/point | net |
  |---|---|---|---|
  | 512 | 239 M | 101.3 | 0.96 GB/s |
  | 1024 | 360 M | 66.5 | 1.44 GB/s |
  | 2048 | 585 M | 40.8 | 2.35 GB/s |
  | 4096 (default) | 689 M | 34.6 | 2.76 GB/s |
  | 16384 | 786 M | 29.8 | 3.14 GB/s |

  Shrinking the pool to fit L3 therefore **backfires**: the only way to shrink it is
  smaller blocks, and the per-message cost on the single service thread grows faster than
  the cache residency helps.
- **`--swc` is nearly irrelevant** above ~128: `--block 1024 --swc 128` and
  `--block 1024 --swc 512` are within noise (66.9 vs 66.5 ns/point).  At np 1 a `--swc`
  sweep is monotone but shallow (32 -> 31.5 ns, 128 -> 24.6, 512 -> 20.9, 1024 -> 20.6);
  per-line overhead, not L1 residency, is what shrinking it costs.
- `--inbox` and `--n-recv` do nothing measurable.

### Ruled out

- **NUMA placement of the pool.**  `router.hpp:558` allocates and `memset`s the whole pool
  from thread 0 (`:628`), so first touch puts all of it on one node.  Real, but not the
  limiter: `numactl --interleave=all` is slightly *worse* (19.2 vs 18.4 ns/point at np 1).
- **Private write-combining lines overflowing L1d** (F x swc x 16 B = 48 KB > 32 KB).  The
  `--swc` sweep goes the other way; see above.

## 3. The benchmark itself distorts the result

`examples/router_bench.cpp:77` is

```c
int d = (int) (x % (u64) F);
```

a **64-bit hardware divide per point**.  `arith.divider_active` = 54 G of 509 G cycles, and
this single line is the top entry in the np 1 sender profile (24 % of sender cycles).
Replacing it with a multiply-shift (`((x >> 32) * (u64) F) >> 32`) in an out-of-tree
variant isolates the terms at np 1, 25 + 6:

| sender inner loop | ns/point |
|---|---|
| PRNG only (`raw`) | 5.6 |
| PRNG + `x % F`, no push | 15.4 |
| PRNG + multiply-shift, no push | 6.6 |
| PRNG + multiply-shift + `Router_Push` | 18.4 |
| PRNG + `x % F` + `Router_Push` (as shipped) | 21.0 |

So the modulo costs **8.8 ns/point on its own**, and `Router_Push` costs **11.8 ns/point**
when measured against a divider-free loop.  In the shipped benchmark the push only *appears*
to cost 5.6 ns because it hides behind the divider's latency.  The two largest items partly
overlap, which means **the shipped "20.9 ns/point" both understates the push path and
charges the router for the benchmark's division.**  Worth fixing in the benchmark: it is
harness cost, and it changes which hotspot you see.

Also: **`--local-only` is dead.**  `examples/router_driver.hpp:34/104` parses it into
`RouterArgs::local_only`, and `router_bench.cpp` never reads it.  The node-local
measurements above come from a patched variant that picks `d` in
`[rank * per_node, (rank + 1) * per_node)`.

## 4. Where the local path's ~5.3 ns/point goes

Marginal cost of `Router_Push` with node-local destinations at np 2 is 11.9 - 6.6 =
**5.3 ns/point**.  At np 1 (both sockets, no MPI) the sender profile splits as: 32 % the
body of `Router_Push`, 7.4 % the 8 KB `memcpy` from the private line into the block, 2.9 %
the `fetch_add` on `dest[]` plus the `n_valid` release store, 2.1 % `cpu_relax`.

Machine-level at np 1: IPC 1.73, `resource_stalls.sb` (store buffer full) 9.3 % of cycles,
`cycle_activity.stalls_mem_any` 21.7 %, `machine_clears.memory_ordering` negligible, and
**42 GB/s of DRAM traffic — 43 bytes of DRAM per 16-byte point** (`uncore_imc` CAS counts),
the block pool's RFO reads plus the receivers' reads.

Note on `Router_Push` (`router.hpp:845-857`): the per-destination count lives in
`line[swc_linesize - 1].val`, i.e. inside the line itself.  Every push therefore touches
**two** cache lines (the slot and the tail line holding the count) and carries a
load -> store -> load chain through that word.  Lines 847 (`u64 n = count`) and 848
(`line[n].key = a`) alone are 7.4 % and 14.7 % of np 1 sender cycles.  Not measured: whether
moving the count into the `Router_thread` handle helps.

Receivers are not the bottleneck in any configuration (senders' `cpu_relax` at np 1 is 2 %).
Their cycles are ~46 % `memmove` (`Router_Pop` copying out of the block) and ~42 %
`murmur128`, which is the benchmark's own work.

## Caveat

np 2 on a single machine is **not** a cluster.  The inter-rank ceiling measured here is a
kernel cross-socket memcpy; on a real multi-node run it would be the NIC instead.  What
does carry over is the structure: one service thread does all of it, and the free-block
stack couples the remote path's back-pressure to every sender, including those pushing to
local receivers (§1).

## Reproducing

```bash
# the two regimes
mpirun -np 1 --bind-to none build/examples/router_bench --senders 25 --receivers 6
mpirun -np 2 --map-by numa --bind-to numa build/examples/router_bench --senders 12 --receivers 3

# senders' spin sites: record, then resolve the pause instructions
perf record -F 999 --all-user -o r.data <binary> ...      # --all-user hides the service thread's kernel work
objdump -d --no-show-raw-insn <binary> | grep pause       # -> three sites
addr2line -i -f -C -e <binary> 0x<site>                   # -> router.hpp:837 / :769 / bench_receiver
perf script -i r.data --tid=<senders> | awk '{print $NF}' | sort | uniq -c | sort -rn

# the service thread is only visible with kernel symbols
perf record -F 999 -o k.data <binary> ...                 # needs kptr_restrict=0
```

---

# Session 2 -- grdix (256 cores, AMD)

Profiling session, 2026-09-07, on **grdix-13** (Grid'5000, Nancy): 2 sockets x AMD EPYC
9754 (Bergamo, Zen 4c), 128 cores / 256 PU per socket, **256 cores / 512 PU total**, 2
NUMA nodes, but **32 L3 domains of 16 MB** (8 cores + their SMT siblings each), 1 TB RAM,
Open MPI 5.0.7.  Rebuilt on the node (`-march=native` binaries do not cross generations)
with `cmake -DCMAKE_BUILD_TYPE=release` (`-O3 -DNDEBUG -march=native -ggdb`); `perf`
needed `perf_event_paranoid = -1` and `kptr_restrict = 0` (`sudo-g5k`).

Everything below is stock `examples/router_bench.cpp` as committed -- **no out-of-tree
variants this time**.  Every number is the median of 3-5 runs of one 3 s round; the range
is given where it matters.

**Do not compare ns/point across the two sessions.**  The EPYC 9754 is a 2.25 GHz part:
the senders' raw PRNG loop costs **10.4 ns/point** here against grvingt's 5.6, so a
single thread is about half as fast.  Only ratios carry over.

## Summary

The router does not scale to this node on its shipped defaults.  There are two separate
ceilings, and **both are reachable with existing knobs**:

| configuration | routed | ns/point |
|---|---|---|
| np 1, 50 senders + 12 receivers | 1.9 G pts/s | 25.7 |
| np 1, 100 + 25 | 2.3 G pts/s | 44.0 |
| np 1, 204 + 51 | 1.3 G pts/s | 179.3 |
| np 1, 204 + 51, `--block 65536` | **2.8 G pts/s** | **84.8** |
| np 2 (one rank per NUMA node), 50 + 12, defaults | 0.28 G pts/s | 345.8 |
| np 2, 50 + 12, `--credit 64 --block 65536` | **2.1 G pts/s** | **44.1** |

Two headlines.  At np 1 throughput **peaks around 100 senders and then falls**: 204
senders route *less* than 50 do.  At np 2 throughput is **pinned at 278 M pts/s no matter
how many senders push** -- and the defaults `--credit 4 --block 4096` are leaving a
**7.8x** factor on the table.

## 1. The np 1 scaling collapse

Stock, `--bind-to none`, sender:receiver held at 4:1, median of 5 runs:

| senders + receivers | routed | ns/point | range |
|---|---|---|---|
| 12 + 3 | 0.50 G | 24.1 | 23.9-24.1 |
| 25 + 6 | 1.0 G | 24.3 | 24.2-24.7 |
| 50 + 12 | 1.9 G | 25.7 | 25.7-26.1 |
| 100 + 25 | 2.3 G | 44.0 | **35.4-69.6** |
| 152 + 38 | 1.2 G | 121.3 | 112.2-143.3 |
| 204 + 51 | 1.3 G | 179.3 | **112.5-216.6** |

Flat to 50 senders (within 1 % run to run), then a cliff.  Note that the variance
explodes at exactly the point the cost does: below the cliff the benchmark is
reproducible to +/- 1 %, above it a run is worth +/- 30 %.

## 2. The hot spot: the install-wait convoy -- and it is *not* starvation

`perf record -F 499 --all-user`, np 1, 204 + 51, 1.38 M samples.  Threads classified by
role from the inlined call chain (`addr2line -i`), so the shares below are of *all*
samples on the machine unless stated:

| site | share of all samples |
|---|---|
| **`router.hpp:834` -- the install-wait's atomic load** | **37.6 %** |
| **`router.hpp:835` -- its test** | **12.1 %** |
| `router.hpp:837` -- its `pause` | 0.6 % |
| `router.hpp:702` -- `free_pop`, after the `blk_link[head]` load | 6.2 % |
| `router.hpp:703` -- `free_pop`'s CAS | 3.2 % |
| `router.hpp:848` / `:851` -- the body of `Router_Push` | 3.1 % / 2.0 % |

The 204 senders are 81.0 % of all samples, so in *steady-state sender cycles* (dropping
the 8.7 % spent in the `raw` phase) that is **~67 % in the install wait** and **~13 % in
the free-stack CAS**.  For contrast, the same measurement at 50 senders puts
`atomic_base.h:501` at **0.21 %** and `router.hpp:835` at **0.06 %** -- the convoy simply
does not exist below the cliff.

**The mechanism is not the one session 1 describes.**  Session 1 §1 attributes the convoy
to the free-block stack running dry because blocks are held by in-flight MPI sends.  Here
there is no MPI at all (np 1, `net 0 msgs/s`), and the stack never empties:

| site | samples | share |
|---|---|---|
| `take_fresh:769` -- the "wait for a block" `pause` | **0** | 0.000 % |
| `free_pop:699` -- the `head == ROUTER_NONE` branch | 79 | 0.006 % |

Zero.  The pool is sized with `slack = 32 * senders` (`router.hpp:551`), 10226 blocks /
671 MB at this size, and it is never exhausted.  What serializes the sealers is the
**CAS on the single `free_top` word** itself: 332 K seals/s, each a contended
compare-exchange plus a dependent `blk_link[head]` load, from 204 threads spread over 32
L3 domains and 2 sockets.  A sealer stuck in that loop leaves its destination's `next` in
the "installing" state (`k > L`), and every other sender pushing to that destination
convoys behind it at `router.hpp:834-837`.

So session 1's *structural* claim survives -- the single global free stack couples every
sender through head-of-line blocking -- but the trigger is **CAS serialization under
thread count**, not back-pressure from the remote path.  At 256 cores it appears with no
network involved.

Receivers are not the bottleneck: 66 % of their cycles sit at `router.hpp:893`, the
`cur_blk == ROUTER_NONE` test guarding `inbox.pop`, i.e. they are spinning on an empty
inbox.  The service thread is 0.4 % of samples (1 thread of 256) and idle at np 1.

## 3. Ruled out

- **Socket / NUMA crossing.**  At a constant 126 threads, one core each: all on socket 0
  is **46.1** ns/point, split 63/63 over both sockets is **46.7**.  Irrelevant.  (Run with
  `OMP_PROC_BIND=close OMP_PLACES=cores` under `taskset -c 0-127` vs `-c 0-62,128-190`.)
- **Free-stack starvation** -- see the zero-sample table above.
- **Receivers** -- starved, spinning on an empty inbox.
- **The service thread** -- at `--block 32768` it handles 8x fewer blocks and throughput
  is unchanged, so its block rate is not what caps np 1.

## 4. Knobs at np 1

`--block` is the lever, because it is what sets the seal rate and therefore how many
threads are in the `free_top` CAS at once (204 + 51, median of 3):

| `--block` | routed | ns/point |
|---|---|---|
| 4096 (default) | 1.3 G | 166.1 |
| 8192 | 2.3 G | 88.9 |
| 16384 | 2.4 G | 84.9 |
| 32768 | 2.6 G | 83.3 |
| 65536 | 2.8 G | 84.8 |

`--dests` (the virtual fan-out, which shortens the convoy without changing the seal rate)
helps less and saturates: 51 -> 162.4, 102 -> 160.3, 204 -> 157.8, 408 -> 125.9,
816 -> 97.2 ns/point.  **The two do not compose**: `--block 32768 --dests 816` is 90.2,
no better than `--block 32768` alone.  Past that point the router is doing real work --
at `--block 65536` the profile is 33 % `memmove` (the line-to-block copy and the
receiver's copy-out) with the convoy down to 21 %.

**Big blocks only pay when contended.**  At 50 + 12, `--block 65536` is *worse* than the
default (28.9 vs 25.7 ns/point): the pool becomes 2.7 GB and goes cold.

## 5. np 2 (one rank per NUMA node): one saturated service thread

`--map-by numa --bind-to numa`, uniform destinations, defaults:

| senders + receivers per rank | routed | ns/point |
|---|---|---|
| 12 + 3 | 289 M | 81.8 |
| 25 + 6 | 277 M | 175.6 |
| 50 + 12 | 280 M | 343.4 |
| 100 + 27 | 277 M | 651.2 |

**Flat.** Throughput is fixed at ~278 M pts/s and the per-point cost is just that ceiling
divided among however many senders you start.  Per rank: 1.11 GB/s of inter-rank traffic,
17.1 K messages/s of 65600 bytes.

The service thread is saturated, and session 1 is right that the reported duty cycle
(here "4 % busy") is a fraction of *turns*, not of cycles.  With kernel symbols, the two
service threads' own profile is **61 % kernel**:

- **54.4 % `_copy_to_iter`**, reached through `process_vm_rw_core` -- Open MPI's `sm` BTL
  using CMA single-copy (this is the AMD equivalent of session 1's `rep_movs_alternative`);
- 11.9 % `ompi_request_default_test_some`, 8.6 % `ompi_request_check_same_instance`,
  1.7 % `PMPI_Testsome`, 1.4 % `opal_progress`;
- the rest per-transfer overhead: `__get_user_pages`, `gup_put_folio`, `mmput`.

It collects 4119 samples against an average of 4420 for the (100 %-spinning) worker
threads, i.e. **~93 % of a core**.  `sys` time for the run is 5.37 s.

### The defaults cost 7.8x

| configuration (np 2, 50 + 12) | routed | ns/point |
|---|---|---|
| defaults (`--credit 4 --block 4096`) | 277 M | 345.8 |
| `--credit 16` | 780 M | 126.9 |
| `--credit 64` | 777 M | 126.3 |
| `--block 16384` | 705 M | 120.2 |
| `--block 65536` | 968 M | 65.4 |
| **`--credit 16 --block 65536`** | **2.1 G** | **44.1** |
| `--credit 64 --block 65536` | 2.1 G | 44.7 |

Unlike at np 1, here the two knobs **do** compose, and together they take np 2 up to the
np 1 ceiling.  In the tuned configuration the inter-rank path carries **8.29 GB/s** per
rank in 7.9 K messages/s of 1048640 bytes.

With the knobs fixed, np 2 saturates gracefully rather than collapsing: 25 + 6 gives
22.7 ns/point, 50 + 12 gives 44.1, 100 + 27 gives 84.0 -- all at ~2.0-2.1 G pts/s.

## 6. Where this contradicts session 1

1. **The convoy's trigger.**  Session 1 §1 blames an empty free stack (back-pressure from
   in-flight MPI sends).  Here the stack is never empty -- literally zero samples in
   `take_fresh`'s wait -- and the convoy is 67 % of sender cycles anyway, at np 1 with no
   MPI.  The cause at scale is CAS serialization on `free_top`, not starvation.
2. **The inter-rank ceiling is a message-rate limit, not a copy-bandwidth limit.**
   Session 1 §2 concludes that "the router's 2.7 GB/s *is* the achievable rate for a
   single thread copying cold, cross-socket 64 KB blocks", from a standalone streaming
   test.  On this node the same single service thread moves **1.11 GB/s at 65600-byte
   messages and 8.29 GB/s at 1048640-byte messages** -- 7.5x from message size alone.
   What is saturated at the default block size is per-message overhead, not the copy.
3. **`--credit 4`.**  Session 1 measures +13 % from raising it.  Here it is worth
   **2.8x** and is the single largest one-flag win at np 2.
4. **The `x % F` divider.**  Session 1 §3 makes it the top np 1 sender entry (24 % of
   sender cycles) and warns it distorts the result.  That is a Skylake property: the one
   `divq` in the binary (`router_bench.cpp:77`, file offset `0x4394`) draws **0.003 %** of
   samples at 204 senders and does not reach the top 15 source lines at 50 senders either.
   On Zen 4c it is free.  **Session 1 §3's warning does not carry over to this machine**
   -- but the modulo is still there, so it will bite again on Intel.
5. **Dead options.**  Session 1 notes `--local-only` is parsed and never read.  So are
   **`--no-bind` and `--cache-level`**: `router_bench.cpp` binds nothing at all
   (`examples/router_driver.hpp:32-34` declare all three and the parser sets them; `router_bench.cpp` reads none).  At 32
   cores that hardly showed.  At 256 it is the main source of run-to-run noise -- at
   100 + 25, one unbound run gave 58.8 ns/point against 37.5 for the same run under
   `OMP_PROC_BIND=spread OMP_PLACES=cores`.  Binding is also what makes the numbers in
   this session reproducible.

## 7. What to change

- **Raise the defaults.**  `credit = 4` and `block_points = 4096` (`router.hpp:170/175`)
  are tuned for a much smaller machine.  `credit >= 16` costs nothing anywhere it was
  measured; `block_points` wants to scale with the team, since big blocks pay only when
  contended (§4).
- **The single `free_top` stack is the structural limit at np 1** (§2).  It serializes
  every sealer on the node through one cache line, and each serialized sealer convoys its
  destination's senders.  Per-thread-group or per-destination free lists, or a per-sender
  cache of a few blocks, would break the coupling; nothing else measured here touches it.
- **Bind the bench's threads**, or delete `--no-bind` / `--cache-level` / `--local-only`
  so the option list stops promising what it does not do.
- The `service N% busy` display is a fraction of turns and reads as idle when the thread
  is at 93 % of a core.  Worth reporting cycles, or renaming.

## Reproducing

```bash
# the np 1 cliff (stock; add OMP_PROC_BIND=close OMP_PLACES=cores to cut the variance)
mpirun -np 1 --bind-to none build/examples/router_bench --senders 50  --receivers 12
mpirun -np 1 --bind-to none build/examples/router_bench --senders 204 --receivers 51

# np 2, defaults vs tuned
mpirun -np 2 --map-by numa --bind-to numa build/examples/router_bench --senders 50 --receivers 12
mpirun -np 2 --map-by numa --bind-to numa build/examples/router_bench --senders 50 --receivers 12 \
       --credit 64 --block 65536

# socket crossing at constant thread count
OMP_PROC_BIND=close OMP_PLACES=cores mpirun -np 1 --bind-to none \
       taskset -c 0-127        build/examples/router_bench --senders 100 --receivers 25
OMP_PROC_BIND=close OMP_PLACES=cores mpirun -np 1 --bind-to none \
       taskset -c 0-62,128-190 build/examples/router_bench --senders 100 --receivers 25

# the spin sites, resolved through the inline chain rather than by hand
objdump -d --no-show-raw-insn build/examples/router_bench | grep pause
addr2line -i -f -C -e build/examples/router_bench 0x<file offset>
perf record -F 499 --all-user -o r.data -- mpirun ...
perf report -i r.data --stdio --no-children -g none --sort srcline

# the service thread is only visible with kernel symbols (kptr_restrict=0)
perf record -F 499 -o k.data -- mpirun -np 2 ...
```

---

# Session 3 -- grvingt, the router after `4a706a1` ("router with thread placement")

Re-run of the campaign on **grvingt-12** again (Grid'5000, Nancy) -- the *same machine as
session 1*, so session 1's numbers are a direct A/B and session 2's (AMD, 256 cores) are
not.  2 sockets x Intel Xeon Gold 6130, 32 cores / 64 PU, 2 NUMA nodes, 2 L3 domains of
22 MB, Open MPI 5.0.7.  `cmake -DCMAKE_BUILD_TYPE=release`, `perf_event_paranoid = -1`,
`kptr_restrict = 0` (`sudo-g5k`).

**Both binaries are stock**, built from the same tree: "OLD" is `a3b5ded` (the code
sessions 1 and 2 measured, in a `git worktree`), "NEW" is `4a706a1`.  Every number is the
median of 3 rounds of 3 s unless stated.  OLD is run `--no-bind`-equivalent (it pins
nothing); NEW pins by default -- §4 shows that makes no difference here.

What changed in `4a706a1`, of the things session 2 §7 asked for: the free stack is now
popped **32 blocks at a time** (`ROUTER_BATCH`, `free_pop_many`) into a **per-sender
cache**, which is exactly the "per-sender cache of a few blocks" recommended there.  The
defaults `credit = 4` and `block_points = 4096` were **not** raised, and the bench still
does not read `--local-only`... it does now, but `x % F` is still there (§6).

## Summary

**The convoy is fixed, and on this machine that buys nothing, because this machine never
had it.**

| claim | session 1/2 | now |
|---|---|---|
| senders in the install wait (`stage_line`) | 28.0 % (S1, np2) / 49.7 % (S2, np1 204 senders) | **0.00 %** |
| senders in the free stack (`take_fresh` / `free_pop`) | 13.5 % (S1) / 9.4 % (S2) | **0.00 %** |
| `--dests`, the convoy-shortening knob | helps 1.7x (S2) | **hurts**, monotonically (§3) |
| lossy vs lossless | -- | **identical** (35.8 vs 35.9): senders never wait |
| np1 throughput vs OLD, 12..54 senders | -- | **3-6 % worse** (§2) |
| np1 throughput vs OLD, 56 senders | -- | ~4 % better |

The fix is real and verifiable, but its target regime -- hundreds of threads convoying on
one `free_top` cache line -- does not exist on 32 cores.  What is left here is its cost:
a bigger, colder block pool (§4), and a **new bistable collapse at small `--block`** that
OLD did not have (§5).

## 1. The convoy is gone

`perf record -F 499 --all-user`, np 1, 50 senders + 12 receivers, 152 K samples, source
lines summed by region of `include/router/workers.hpp`:

| region | share of all samples |
|---|---|
| free stack (`free_push` / `free_pop_many` / `refill`, `:13-80`) | **0.00 %** |
| the install wait (`stage_line`'s `k > L` spin, `cpu_relax` at `:153`) | **0.00 %** |
| `Router_Push` body (`:157-175`) | 23.07 % |
| everything else in `workers.hpp` | 0.60 % |

The only free-stack samples anywhere are `workers.hpp:18/19` at 0.00 %.  For comparison
session 1 measured 28.0 % + 13.5 % and session 2 measured 49.7 % + 9.4 % in exactly these
two regions.

Two independent corroborations, no profiler needed:

- **`--dests` has flipped sign.**  Session 2 §4 used it to shorten the convoy and got
  1.7x.  Here it only adds destinations to stage into: 12 -> 35.5, 24 -> 36.8, 48 -> 38.9,
  96 -> 40.2, 192 -> 42.7 ns/point.
- **Lossy is not faster than lossless** (35.8 vs 35.9 ns/point, 0.76 % dropped).  If
  senders were waiting, dropping instead would show.

The rest of the profile is real work: `router_bench.cpp:79` (the `x % F` divider, §6) at
23.3 %, `Router_Push`'s two stores into the private line (`workers.hpp:163/164`) at
10.7 % + 5.1 %, then the PRNG, `murmur128` and `memmove`.

## 2. A/B against the old router, np 1

Median of 3, and OLD reproduces session 1's 25+6 number (20.9) exactly, which validates
the harness:

| senders + receivers | OLD ns/pt | NEW ns/pt | NEW routed |
|---|---|---|---|
| 12 + 3 | 17.9 | 19.1 | 0.63 G |
| 25 + 6 | **20.9** (S1: 20.9) | 22.3 | 1.1 G |
| 36 + 9 | 28.2 | 29.1 | 1.3 G |
| 44 + 11 | 32.6 | 33.6 | 1.4 G |
| 50 + 12 | 34.8 | 35.7 | 1.4 G |
| 54 + 9 | 35.3-35.7 | 36.4-36.7 | 1.5 G |
| 56 + 7 | 44.6-46.1 | 43.2-44.8 | 1.3 G |
| 58 + 5 | 60.6-72.9 | 63.7-66.7 | 0.9 G |

Run-to-run spread is under 1 % up to 54 senders for both, so the 3-6 % gap is real.  NEW
wins only at 56 + 7, by ~4 % -- the first pass suggested 20 % there, but repeats put it at
4 %; do not quote the single run.  At 58 + 5 both collapse (too few receivers) and the
difference is noise.

There is **no session-2-style cliff on this node**: throughput rises to ~1.4-1.5 G pts/s
and stays there.  32 physical cores cannot reach the ~100-sender regime where session 2
saw 2.3 G fall to 1.3 G.

## 3. Where the 3-6 % goes: a colder pool

The per-sender cache is charged to the block pool.  Same run, 50 + 12, from the banners:

| | blocks | pool |
|---|---|---|
| OLD | 2570 | 168.6 MB |
| NEW | 4170 | 273.6 MB |

The difference is exactly 50 x 32 = 1600 blocks: `connect.hpp:220` adds
`S * (ROUTER_BATCH + 1)` on top of the old `slack = max(F, 32 * S)`.  Both pools are far
past the 44 MB of L3, so the extra 105 MB is pure cold DRAM, and it shows on the memory
controllers (`uncore_imc` CAS counts, 3 s, 50 + 12):

| | DRAM read + write | per point |
|---|---|---|
| OLD | 149.9 GiB | ~35 bytes |
| NEW | 166.9 GiB | ~42 bytes |

+20 % of DRAM bytes per point for 4 % fewer points.  (Session 1 measured 43 bytes/point at
25 + 6; the shape is the same.)

## 4. Placement is a no-op on this node

`4a706a1` gives the Router its own pinning and grouping.  At 50 + 12 it is worth nothing
measurable: NEW pinned 35.8-36.0, NEW `--no-bind` 35.8-36.0 ns/point.  Expected -- two L3
domains and a flat mask; session 2 §6.5 saw binding matter at 256 cores, not here.
`--group` (4 / 8 / 16 cores per group) is equally flat: 36.1 / 35.9 / 36.0.

The placer's refusal is correct and its message is good: `--senders 50 --receivers 20`
(71 threads for 64 CPUs) is rejected with the three ways out rather than silently
oversubscribing.

## 5. New: a bistable collapse at small `--block`

**This is a regression, and it is the one thing here worth acting on.**  At
`--block 1024`, 50 + 12, NEW lands in one of two modes; OLD is stable:

| | runs |
|---|---|
| NEW | 38.6, 41.3, 41.4, **123.1**, 53.6, **130.3** ns/point |
| OLD | 48.5, 49.0, 52.3 ns/point |

The good mode (~40) beats OLD; the bad mode (~110-130, throughput 1.4 G -> 0.49 G) is 2.5x
worse than OLD.  It is stochastic and it **latches**: 2 s rounds never collapsed in 18
tries, 3-4 s rounds do.  `--block 2048` is stable (36.2-37.0) and the **default 4096 is
rock stable** at every sender count tried -- 25 / 36 / 50 senders, 6 repeats each, all
within 1 %.  So no shipped configuration is exposed today.

The mechanism, from a profile of a caught bad run (108 ns/point), addresses resolved
through the inline chain with `addr2line -i`:

| site | samples (of 152 K) |
|---|---|
| **`refill`, `workers.hpp:78` -- the retry `free_pop_many` in the "wait for a block" loop** | **49.9 K (33 %)** |
| **`free_push`, `workers.hpp:19` -- the CAS, reached from `Router_Pop`** | 9.8 K (6.4 %) |
| `router_bench.cpp:79` (the divider) | 11.8 K |

Region sums for that run: free stack / refill **4.78 %** (against 0.02 % in the good mode),
install wait 0.01 %, `Router_Push` down to 6.43 %.  So:

1. 50 senders x 32 cached blocks = 1600 blocks are **parked in sender caches**, out of
   4170.  The pool is sized to cover that, but the *distribution* is not: a sender holding
   32 idle blocks starves one that has none.
2. The free stack empties.  Starved senders spin in `refill`'s `while (n == 0)` loop --
   and that loop calls `free_pop_many`, i.e. it **re-loads the contended `free_top` on
   every iteration** rather than backing off on a private word.
3. That load traffic on `free_top` is what the receivers' `free_push` CAS has to win to
   *return* blocks.  It cannot, so the stack stays empty.  **Positive feedback -- which is
   why the state latches instead of passing.**

Note this is the *same* head-of-line coupling through one global free stack that sessions
1 and 2 both blamed; batching moved it from the steady state into a rarer but far worse
metastable state.  A backoff on a private word in `refill`'s wait, or reserving blocks the
receivers can always push into, would break step 3.

## 6. Unchanged from session 1

- **`x % F` is still in the bench** (`router_bench.cpp:79`, now covering `--local-only`
  too) and is still the top single source line at np 1: **23.3 %** of samples.  Session 1
  §3's warning stands on Intel, and every np 1 ns/point here still includes it -- the
  `raw` line puts the sender's PRNG + divider loop at 5.8 ns/point at 25 senders and
  9.5 at 50.  Against that floor NEW's marginal push cost at 50 + 12 is 26.4 ns/point.
- **`--local-only` now works** (it was dead in session 1): np 2, 12 + 3 node-local gives
  20.2 ns/point against 34.9 uniform, consistent with session 1's patched 11.9 + the
  ~8.8 ns divider.

## 7. np 2, for the record

Not the target topology any more -- with the placer forming groups per cache domain, one
rank per node is the intended shape -- but the numbers were taken before that was settled,
and they confirm session 2 §6.3 on Intel:

| np 2, 12 + 3, uniform | ns/point |
|---|---|
| defaults (S1: 34.7) | 34.8 |
| `--credit 16` | 30.5 |
| `--credit 64` | 30.0 |
| `--block 65536` | 34.9 |
| **`--credit 16 --block 65536`** | **23.0** |

The inter-rank path is **unchanged at defaults** (34.8 vs session 1's 34.7, 2.77 GB/s,
42 K msgs/s) -- as expected, since it is the service thread, not the senders, that is
saturated there, and nothing in `4a706a1` touches it.  Session 2's "the knobs compose"
holds here too: 34.8 -> 23.0, a 1.5x that costs one flag pair.

## 8. Everything still works

`ctest` passes all three (`router_np1/2/4`).  The engine, whose `placement.hpp` this commit
rewrote on top of the new shared `router/topology.hpp`, passes all three smoke tests of
CLAUDE.md: PCS (`--n 20`), direct (`--engine direct`), and the `--dp-len-bits 1` resolver
torture test, each finding its golden pair.

## What to change

- **Fix `refill`'s wait loop** (§5): back off on a private word instead of re-loading
  `free_top` every iteration, so a starved sender cannot block the receivers from
  refilling the stack.  This is the only defect found.
- **`ROUTER_BATCH = 32` is a lot to park per sender**, and it is charged twice to the pool
  (`slack` is already `32 * S`, and `connect.hpp:220` then adds `S * (ROUTER_BATCH + 1)`).
  Worth checking whether 8 keeps the session-2 win at 256 cores while giving back the
  3-6 % and the 105 MB here.
- Session 2 §7's other two items are still open: `credit = 4` is still the default and is
  still worth 1.5x at np 2 (§7), and the bench still divides (§6).

## Reproducing

```bash
git worktree add /tmp/old-router a3b5ded && cmake -S /tmp/old-router -B /tmp/old-router/build \
    -DCMAKE_BUILD_TYPE=release && make -C /tmp/old-router/build router_bench

# the A/B
mpirun -np 1 --bind-to none /tmp/old-router/build/examples/router_bench --senders 50 --receivers 12 --no-bind
mpirun -np 1 --bind-to none build/examples/router_bench --senders 50 --receivers 12

# the bistable collapse: repeat, 3 s rounds or longer, and watch for ~110 ns/point
for i in $(seq 6); do mpirun -np 1 --bind-to none build/examples/router_bench \
    --senders 50 --receivers 12 --block 1024 --seconds 3 --rounds 1; done

# DRAM per point
perf stat -a -e uncore_imc/cas_count_read/,uncore_imc/cas_count_write/ -- mpirun ...

# the inline chain (PIE: file_off = ip - mmap_vaddr + mmap_pgoff, from --show-mmap-events)
perf record -F 499 --all-user -o r.data -- mpirun ...
perf report -i r.data --stdio --no-children -g none --sort srcline
addr2line -i -f -C -e build/examples/router_bench 0x<file offset>
```

---

# Session 4 -- grdix again, the router after `4a706a1`, at 255 workers

Profiling session, 2026-09-08, on **grdix-8** (Grid'5000, Nancy): the same machine class as
session 2 -- 2 sockets x AMD EPYC 9754 (Bergamo, Zen 4c), 256 cores / 512 PU, 2 NUMA
nodes, **32 L3 domains of 16 MB** (8 cores each), 1 TB RAM, Open MPI 5.0.7.  Stock
`examples/router_bench.cpp` at `56f4fb8`, `cmake -DCMAKE_BUILD_TYPE=release`; `perf`
needed `perf_event_paranoid = -1` and `kptr_restrict = 0` (`sudo-g5k`).

Session 2 measured the router *before* the per-sender block cache; session 3 measured the
cache on 32 cores, where its target regime does not exist.  This is the cache at 256
cores, which is what it was written for.  The starting point was two configurations, both
**255 workers + 1 service = 256 threads on 256 cores**, differing only in the split:

| configuration | routed | ns/point | dropped |
|---|---|---|---|
| `--senders 150 --receivers 105 --lossy` | 3.1-3.2 G pts/s | 48 | 0.02 % |
| `--senders 192 --receivers 63 --lossy` | **1.3 G** or **2.7 G** | 88 or 82 | **40 %** or 0.01 % |

The second one is bistable.  Explaining that, and the ceiling both of them run into, is
what this session is about.

## Summary

Three findings, in decreasing order of how much they cost:

1. **The node's ceiling is a *block* rate, not a point rate**, and it is the single
   `free_top` word -- specifically the senders' **pop** side (§3).  Across the whole
   campaign -- senders 85..234, receivers 21..170, `--block` 1024..4096 -- the router
   handles **460-830 K blocks/s** and nothing else.  At the default 4096-point block that is
   2-3 G pts/s, and it is why `--block` is the only lever that matters: raising it divides
   the free-stack traffic per point.  Senders spend **36 % of their cycles** in
   `free_pop_many` in the healthy state and **74 %** where the ceiling binds.  Chaining the
   receivers' frees, which a microbenchmark said was worth 17x, was implemented and is
   worth **zero**: the contention is on the pop side.
2. **Lossy mode has a latching collapse that lossless does not** (§2).  Drops return blocks
   to the free stack directly, so the senders never wait, the receivers' inboxes stay full,
   the write-to-read delay grows past the lifetime of a line in cache, and the receivers'
   copy-out slows by 2x -- which doubles the delay.  Positive feedback.  A 3:1 split sits
   right on the tipping point; 1.4:1 has enough receiver headroom to stay out of it.
3. **The placer stacked both round-robin passes on the low domains** (§5), doubling up
   cores at one end of the machine while leaving cores idle at the other.  Fixed here in
   three lines; worth **+8 %** at np 1 and **+15 %** at np 2.

Session 3's fix works and is visible: the install-wait convoy that was 49.7 % of sender
cycles in session 2 is **gone** (0.00 % here too), and `take_fresh`'s starvation wait with
it.  What is left is the free stack it pops from.

**Do not compare ns/point with session 3** (grvingt, 32 Intel cores).  The senders' raw
PRNG loop is 10.4 ns/point here, as in session 2.

## 1. The split, swept at a constant 255 workers

`--lossy`, `--seconds 2 --rounds 3`, the last round of each:

| S + R | routed | ns/point | dropped | per receiver |
|---|---|---|---|---|
| 234 + 21 | 471 M | 84.0 | 83.8 % | 22 M/s |
| 223 + 32 | 715 M | 81.5 | 73.9 % | 22 M/s |
| 213 + 42 | 916 M | 88.8 | 63.9 % | 22 M/s |
| 204 + 51 | 1.1 G | 93.0 | 53.5 % | 22 M/s |
| 192 + 63 | 2.6 G | 81.8 | 0.01 % | 41 M/s |
| 171 + 84 | 1.8 G | 97.8 | 7.6 % | 21 M/s |
| 150 + 105 | 3.2 G | 47.9 | 0.02 % | 30 M/s |
| **127 + 128** | **3.4 G** | **37.1** | 0.01 % | 27 M/s |
| 105 + 150 | 3.0 G | 34.6 | 0.01 % | 20 M/s |
| 85 + 170 | 3.2 G | 27.0 | 0.00 % | 19 M/s |

Two populations, not a curve.  Every configuration is either **healthy** (drops in the
1e-4 range, 3-3.4 G pts/s) or **collapsed** (half to five sixths of the points thrown
away, 0.5-1.8 G), and the giveaway is the last column: a receiver delivers **41, 30, 27,
20 M pts/s** when healthy and **21-22 M pts/s** in every collapsed run, whatever the
sender count.  The collapsed state has its own, lower, per-receiver rate.

`192 + 63` belongs to neither: it is the boundary, and it lands on one side or the other
per round.  Five rounds of 3 s, twice over, gave 1.4/1.3/1.3/1.3/1.3 G both times; at 2 s
and at 4 s the first round collapsed and the rest recovered to 2.6-2.7 G.  **Round 0
almost always collapses** (the senders start with 32 cached blocks each and the receivers
with nothing, so the queues fill before the receivers get going), and whether it climbs
back out is a coin toss.

## 2. The collapse: lossy mode removes the only flow control there is

Every drop in the collapsed runs is `DROPPED_RECV`, which `service.hpp:96` raises when
`receivers[local]->inbox.push()` fails -- the receiver's inbox is **full**.  So the
receivers are the constraint.  They are not, however, idle: `perf record -a -F 299`,
threads separated by CPU (the placer pins them, so `perf report -C` is an exact role
filter):

| | healthy (lossless) | collapsed (lossy) |
|---|---|---|
| receiver cycles in `memmove` (`Router_Pop`'s copy-out) | ~50 % | **79 %** |
| receiver cycles in the pop loop / empty inbox (`workers.hpp:207/209`) | 22.6 % | **0 %** |
| receiver IPC | 0.16 | **0.04** |
| points delivered per receiver | 38 M/s | 21 M/s |

The collapsed receiver never waits and never spins: it is inside `rep movsb` at IPC 0.04.
Same code, same instruction stream, same bytes, **2x the cycles per byte**.  That is
memory, and the counters say which memory -- L2 prefetches that miss both L2 and the
local L3, per point delivered:

| regime | per receiver | `l2_pf_miss_l2_l3` per point |
|---|---|---|
| lossless, `--block 4096` | 38 M/s | 0.103 |
| lossy collapsed, `--block 4096` | 21 M/s | **0.174** |
| lossy, `--block 16384` | 55 M/s | 0.100 |

**1.7x the far traffic per point, for 1.8x the cost.**  The block a collapsed receiver
reads has gone cold: by the time it is read, the lines the sender wrote are no longer in
the sender's 16 MB L3 domain.

Why lossy and not lossless.  A dropped block is `release`d immediately (`service.hpp:104`),
straight back through the service's stash to the free stack, so **a sender never runs out
of blocks and never waits**: it pushes at full rate no matter how far behind the receivers
are.  The inboxes therefore sit full, and the queueing delay between the write and the read
is `in-flight bytes / drain rate`.  Slow the drain and the delay grows, which cools the
lines further, which slows the drain: the state latches instead of passing.  In lossless
mode there is nothing to break the coupling -- a sender out of blocks waits in `refill`, so
the senders can never outrun the receivers, the queues stay short, and the collapse is
**unreachable**:

| 192 + 63, 4 rounds of 3 s | routed |
|---|---|
| lossless | 2.3 / 2.4 / 2.4 / 2.4 G, 0.000 % dropped, every time |
| lossy | 1.3-1.4 G latched, or 2.6-2.7 G, per round |

The queue depth is the control parameter, and `--inbox` sets it (192 + 63 lossy, 3 rounds
of 3 s):

| `--inbox` | rounds |
|---|---|
| 4 | 2.3 G, 2.6 G, 2.7 G (1.2 % dropped) |
| 8 | 2.3 G, 2.6 G, 2.6 G (0.06 %) |
| 16 | 2.2 G, 2.6 G, 2.6 G (0.04 %) |
| 32 | 2.2 G, 2.6 G, 2.6 G (0.02 %) |
| **64 (default)** | 1.4 G, 2.8 G, 2.8 G -- **bistable** |
| 128 | 1.4 G, 1.5 G, 1.5 G (42 % dropped) -- latched |
| 256 | 1.4 G, 1.5 G, 2.6 G |

Shallower inboxes do not fix it (`--inbox 32` still collapsed in a 1-round run), they only
make the latch harder to reach.  What actually removes it is either lossless mode or enough
receivers.

**Ruled out for the collapse:** NUMA.  The pool is 1144 MB and **87 % of it is on NUMA
node 0** (`/proc/PID/numa_maps` during a run -- `router_alloc` memsets from thread 0, so
first touch puts it all there), and `numactl --interleave=all` is worth a real **+17 %**
lossless (2.4 -> 2.8-2.9 G) but does not stop the lossy collapse (2.4 G, 1.7 G, 1.3 G,
1.3 G).  Worth fixing on its own account; not the trigger.

## 3. The ceiling: one `free_top` word, and it is a block rate

The healthy configurations all stop in the same place, and it is not where the point rate
suggests.  Sweeping `--block` while reading the block rate off `pushed / block_points`:

| | `--block` 1024 | 2048 | 4096 | 8192 | 16384 | 32768 |
|---|---|---|---|---|---|---|
| 192 + 63, routed | 650 M | 1.2 G | 1.3 G | 1.9 G | 4.1 G | 4.4 G |
| 192 + 63, **blocks/s** | 635 K | 586 K | 537 K | 390 K | 311 K | 174 K |
| 150 + 105, routed | 705 M | 1.3 G | 2.8 G | 2.5 G | 2.6 G | 2.7 G |
| 150 + 105, **blocks/s** | 688 K | 635 K | 684 K | 305 K | 159 K | 82 K |

Below 4096 points per block the point rate falls exactly as the block size does, at a
**pinned 590-690 K blocks/s**.  Over the whole session -- 10 splits, 6 block sizes -- the
block rate never leaves 460-830 K/s while the point rate moves by a factor of seven.  That
is the signature of a per-block serialization point.

It is the free stack.  The lossless sender profile, by source line
(`include/router/workers.hpp`):

| site | share of sender cycles |
|---|---|
| `:52-55` -- `free_pop_many`'s walk of the 32 links | **30.5 %** |
| `atomic_base.h:501` -- its `compare_exchange_weak` | 5.4 % |
| `:164/:167` -- the body of `Router_Push` | 6.7 % |
| `tools.hpp:133-137` -- the PRNG | ~17 % |
| the install wait (`stage_line`, `k > L`) | **0.00 %** |
| `take_fresh`-style starvation wait | **0.00 %** |

Session 3's fix holds at 256 cores: the convoy and the starvation wait are gone.  What the
per-sender cache did was move the cost from *waiting on* the stack to *fighting for* it --
`free_pop_many` loads `free_top`, walks up to 32 `blk_link` entries (a dependent load
chain over lines other threads are writing), then CASes, and on failure walks again.

Extracted and run on its own -- the same word, the same `blk_link` array, threads popping
a batch of 32 and pushing the 32 back one at a time, which is exactly the router's mix
(one batch-pop per 32 blocks from a sender, one single push per block from a receiver):

| threads | blocks/s through `free_top` | fairness spread |
|---|---|---|
| 1 | 147 M | 1x |
| 8 | 4.5 M | 1x |
| 32 | 1.2 M | 5x |
| 64 | 722 K | 12x |
| 192 | 321 K | 44x |
| **255** | **254 K** | **54x** |

The saturated ceiling is the same order as the router's observed 460-830 K blocks/s, and
the 54x spread between the luckiest and unluckiest thread is visible in the bench too: at
`--block 1024` it reports 2842 ns/point against an aggregate of 311, which is what a mean
of `1/n_i` looks like when some senders are starved.

Both candidate fixes, measured in that same harness at 255 threads:

| | blocks/s | vs today |
|---|---|---|
| today: one stack, one CAS per block returned | 229 K | -- |
| chain the receiver's frees (one CAS per 32) | 3.84 M | 17x |
| **shard per L3 domain (32 stacks)** | **75 M** | **327x** |
| both | 237 M | 1000x |

### Which side of the stack, and what the harness gets wrong

**It is the pop side, and the harness's 17x for chaining is a mirage.**  Chaining the
receiver's frees was implemented in the router (a per-receiver batch of `ROUTER_BATCH`
blocks, `free_push_chain`ed when it fills and whenever the receiver finds nothing left to
read, so it never sits on a block a lossless sender waits for; `router_test` passes,
including with `ROUTER_PARANOID`'s block audit extended to count what a receiver holds) and
it is worth **nothing at all**:

| 192 + 63 | shipped | chained frees |
|---|---|---|
| lossy, `--block 1024` | 623 / 648 / 640 M | 611 / 639 / 645 M |
| lossy, `--block 2048` | 1.2 / 1.3 / 1.3 G | 1.2 / 1.2 / 1.2 G |
| lossless | 2.3 / 2.8 / 2.8 G | 2.6 / 2.7 / 2.7 G |
| 150 + 105, lossy | 2.7 / 3.0 / 3.0 G | 2.6 / 3.0 / 3.0 G |

The reason is in the profile of the regime where the ceiling actually binds -- `--block
1024`, 574 M pts/s, 564 K blocks/s, the service idle at 2.3 M turns/s:

| role | site | share of that role's cycles |
|---|---|---|
| **senders** | `workers.hpp:52-55` -- `free_pop_many`'s walk | **59.5 %** |
| **senders** | `atomic_base.h:501` -- its `compare_exchange_weak` | **14.8 %** |
| senders | `tools.hpp:133-137` -- the PRNG | 9.5 % |
| receivers | `workers.hpp:207/209` -- waiting on an empty inbox | 51.2 % |
| service | idle (`MPI_Testsome` on nothing, 0 % of turns busy) | -- |

**Three quarters of the senders' cycles are in `free_pop_many`.**  The receivers are
starved and the service is asleep: the whole node is 192 senders queueing for one word.
The push side that chaining fixes is 63 threads doing one CAS per block, and it was never
the constraint.

Which of the walk and the CAS?  `ROUTER_BATCH` sets both -- the walk is that many
dependent `blk_link` loads, the CAS is amortized over that many blocks -- and sweeping it
separates them:

| `ROUTER_BATCH` | `--block 1024` | `--block 4096` |
|---|---|---|
| 1 | 333 M | 1.2 G, 0.007 % dropped |
| 4 | 330 M | 1.5 G, 0.017 % |
| 8 | 473 M | **2.0 G, 0.008 %** |
| **32 (today)** | 611 M | 1.4 G, **38 % dropped** |
| 64 | 677 M | 1.4 G, **40 % dropped** |

Monotone in the batch: a **longer** walk is better, so the walk is not what costs -- the
CAS on the one word is, and the only thing that helps is doing fewer of them.  Sharding
does exactly that, by cutting the contenders per word from 192 to the six senders of a
cache domain.

The right-hand column is the sting in the tail, and §2 explains it: **the free stack is
currently also the flow control.**  At `ROUTER_BATCH` 8 the senders are throttled enough
that they cannot outrun the receivers, and `192 + 63` is a stable 2.0 G with nothing
dropped; at 32 and 64 they can, and it collapses to 1.4 G with two points in five thrown
away.  Relieving the stack without giving lossy mode some other back-pressure trades a
throughput ceiling for wasted work.

And in the router itself, relieving the stack by raising the block size takes `192 + 63`
from 1.3 G to **4.4-4.6 G pts/s at 35-36 ns/point** (still 22-27 % dropped: with the stack
out of the way the receivers become the limit).  Note the sign flips with the split, as in
session 2 §4 -- at `127 + 128`, where the stack is not the binding constraint, big blocks
*hurt*: 3.3 G at 4096 against 2.6 G at 16384.

**Ruled out for the ceiling:**

- **The service thread.**  38 % of its cycles at np 1 are Open MPI bookkeeping
  (`ompi_request_check_same_instance` 25 %, `ompi_request_default_test_some` 12.6 %) --
  `MPI_Testsome` over 32 always-posted receives on a run with **no peers and zero
  messages**.  Dead work, but not the limit: `--n-recv 1` against `--n-recv 32` at
  `--block 1024` gives 606-648 M pts/s either way.  Worth skipping when `n_nodes == 1`
  regardless.
- **Receiver capacity per se.**  One receiver alone absorbs **592 M pts/s** (64 senders,
  1 receiver, 82 % dropped).  Per-receiver throughput falls as receivers are added
  (384, 300, 250, 144, 87, 49 M/s at R = 2, 4, 8, 16, 32, 63): a shared limit, not a
  per-thread one.
- **The inbox ring.**  `ring.hpp` already keeps head and tail on separate cache lines.

## 4. Where the time goes, in one place

`192 + 63`, lossless, the healthy state, all 256 threads:

| role | threads | what they are doing |
|---|---|---|
| senders | 192 | 36 % free stack (`free_pop_many` walk + CAS), 17 % PRNG, ~7 % the push body, the rest `memmove` into blocks |
| receivers | 63 | ~50 % `memmove` out of blocks, 23 % waiting on an empty inbox, ~5 % `murmur128` |
| service | 1 | 38 % `MPI_Testsome` on an idle network, the rest `sweep`/`place`/the stash |

Collapsed (lossy), the same team: senders lose the free stack from their profile (they
never starve, drops keep the stack fed) and go ~38 % `memmove`; receivers go **79 %
`memmove` at IPC 0.04** and stop waiting entirely; the service goes to 90 % of its turns
busy, dropping four blocks in ten.

## 5. The placer stacked its two passes (fixed)

`place_pinned` put the service on a core of domain 0 and then started **both** round-robin
passes at domain 0 again, so both remainders landed on the low-numbered domains, on top of
the service:

| | NUMA 0 / NUMA 1 | distinct cores | doubled up | idle |
|---|---|---|---|---|
| 150 + 105, before | **138 / 118** | 246 of 256 | 10 cores | 10 cores |
| 192 + 63, before | 129 / 127 | 255 of 256 | 1 core | 1 core |
| 50 + 12, before | 45 / 18 | 63 | -- | -- |
| **either, after** | **128 / 128** | **256** | **none** | **none** |

At `150 + 105` that is ten cores running two spinning threads each while ten cores at the
other end of the machine sit idle -- and the service thread is one of the pairs, sharing a
core with a sender.  The fix is one cursor for the whole team, started after the service's
domain and carried from the receiver pass into the sender pass, which lands both
configurations on exactly 8 threads per domain (`include/router/placement.hpp`, 3 lines).
Measured, median of the steady rounds:

| configuration | before | after |
|---|---|---|
| np 1, 150 + 105, lossless | 2.9 G, 52.8 ns | **3.1-3.2 G, 47.4-48.2 ns** (+8 %) |
| np 1, 192 + 63, lossless | 2.4 G, 81.0 ns | **2.6-2.7 G, 72.7 ns** (+8 %) |
| np 2, 96 + 31 per rank | 242 M, 774 ns | **280 M, 599 ns** (**+15 %**) |
| np 2, 75 + 52 per rank | 242 M, 555 ns | **279 M, 502 ns** (+15 %) |
| np 1, 50 + 12, lossless | 2.1 G, 23.8 ns | 2.0-2.1 G, 25.5 ns (**-2 to -5 %**) |

The last row is the trade-off and it is worth knowing: for a team far smaller than the
machine, the old placer's accident of packing 45 threads onto socket 0 and 18 onto socket 1
beat an even 32/31 split, because it kept the traffic on one socket.  On the topology the
Router is actually for -- one rank per NUMA node, which it warns about when you do not --
that cannot happen, and there the fix is a clean 15 %.

`ctest` passes all three (`router_np1/2/4`), and the engine passes all three smoke tests of
CLAUDE.md (PCS, `--engine direct`, and the `--dp-len-bits 1` resolver torture test).

## 6. np 2 still ships the wrong defaults

Session 2 §7 asked for `credit` and `block_points` to be raised; they were not, and on this
machine the bill has grown:

| np 2, 96 + 31 per rank | routed | ns/point |
|---|---|---|
| defaults (`--credit 4 --block 4096`) | 280 M | 599 |
| `--credit 16 --block 65536` | **1.9 G** | **81** |

**6.8x, from two flags.**  Unchanged in substance from session 2 §5.

## What to change

- **Shard the free stack per cache domain** (§3).  It is *the* ceiling at np 1 -- 74 % of the
  senders' cycles at the ceiling -- and the router already has the concept it needs: a group
  is senders and receivers in one cache domain, so one free stack per group cuts the
  contenders per word from 192 to 6.  Worth 327x in isolation.  Two things to design: a
  steal or rebalance path (a block moves from the sender's group to the receiver's, so
  shards drift), and where the service's stash spills.  While rewriting it, consider
  dropping the linked list for **an array of block ids per shard with an atomic index**: a
  `fetch_add` is wait-free, needs no walk and cannot fail, and the walk exists only because
  the free list is threaded through `blk_link`.
- **Do NOT bother chaining the receiver's frees** (§3): implemented and measured, worth
  nothing.  The harness that predicted 17x measures a symmetric push/pop load; the router's
  load is 192 senders on the pop side against 63 receivers on the push side.
- **Pair it with real back-pressure for lossy mode** (§2, §3): the free stack is what
  currently keeps the senders from outrunning the receivers, and `ROUTER_BATCH` 8 vs 32 is
  the demonstration -- 2.0 G with nothing dropped against 1.4 G with 38 % dropped.  Fixing
  the ceiling without replacing that throttle just converts it into waste.  `ROUTER_BATCH
  = 8` is a one-constant mitigation available today.
- **Raise the defaults**, again (§6): `credit = 4` costs 6.8x at np 2 with `block_points`.
- **Give lossy mode back some back-pressure** (§2).  Dropping a block returns it to the free
  stack, which is precisely what lets the senders outrun the receivers and cool the pipe.
  A sender whose destination is dropping could be made to slow down, or the queue depth
  bounded in bytes rather than in blocks; `--inbox 16` is a cheap partial mitigation today.
- **First-touch the pool in parallel** (§2): 87 % of 1144 MB lands on one NUMA node, and
  interleaving it is worth 17 %.
- **Skip the MPI progress path when `n_nodes == 1`** (§3): 38 % of the service thread's
  cycles are `MPI_Testsome` over an idle network.  Not a bottleneck, but it is free.

## Reproducing

```bash
# the two configurations, and the bistability (round 0 nearly always collapses)
mpirun -np 1 --bind-to none build/examples/router_bench --senders 150 --receivers 105 --lossy
mpirun -np 1 --bind-to none build/examples/router_bench --senders 192 --receivers 63 --lossy --rounds 5
mpirun -np 1 --bind-to none build/examples/router_bench --senders 192 --receivers 63          # lossless: never collapses

# the block rate is the ceiling: read pushed/block_points, not the point rate
for b in 1024 2048 4096 8192 16384 32768; do mpirun -np 1 --bind-to none \
    build/examples/router_bench --senders 192 --receivers 63 --lossy --block $b; done

# roles separated by CPU, which the placer pins (see the CPUS lines of a placement probe)
perf record -a -F 299 -o r.data -- mpirun -np 1 --bind-to none build/examples/router_bench ...
perf report -i r.data --stdio --no-children -g none --sort srcline -C <sender cpus>

# the receiver's far traffic per point
perf stat -a -C <receiver cpus> -e l2_pf_miss_l2_l3.all,l2_pf_miss_l2_hit_l3.all,cycles,instructions -- mpirun ...

# where the pool's pages are
awk '/N0=/ {print}' /proc/$(pgrep -f '^build/examples/router_bench')/numa_maps

# the free stack on its own (the harness is in the session's scratch, ~70 lines):
# 255 threads, pop a batch of 32 and push the 32 back one at a time, vs chained pushes and vs 32 shards
# -- but see §3: it overstates chaining, because the router's load is asymmetric (192 poppers, 63 pushers)

# which side of the stack costs: sweep the batch, which sets both the walk length and the CAS rate
sed -i 's/ROUTER_BATCH = 32/ROUTER_BATCH = 8/' include/router/common.hpp   # then rebuild
```

---

# Session 5 -- both nodes at `739b7ea`: the free ring works, and how far the memory wall is

Session of 2026-09-08 on **grvingt-11** (2 x Xeon Gold 6130, 32 cores / 64 PU, 2 NUMA nodes,
2 L3 of 22 MB, 187 GB) and **grdix-5** (2 x EPYC 9754 Bergamo, 256 cores / 512 PU, 2 NUMA
nodes, 32 L3 of 16 MB, 1 TB), Grid'5000 Nancy, Open MPI 5.0.7.  Stock tree at `739b7ea`,
`cmake -DCMAKE_BUILD_TYPE=release` in a worktree per cluster; `perf` with
`perf_event_paranoid = -1` and `kptr_restrict = 0` (`sudo-g5k`).  Median of 3 rounds of 3 s.

**np 1 only.**  One rank now holds the whole machine, so the np 2 columns of sessions 1-4 are
not re-measured here.

What changed since session 4 measured this code: `3747c82` made the free stack a bounded MPMC
ring of block ids, `4dd6316` first-touches the pool in parallel, `e75ff5f` gave the receiver a
zero-copy path.  Sessions 3 and 4 asked for the first two by name.

## Summary

1. **Session 4's ceiling is gone.**  `free_pop_many` was 30-59 % of sender cycles there and
   the node was pinned at 460-830 K blocks/s.  The free ring is now **0.13 % (grvingt) and
   0.38 % (grdix) of all samples**, and what is left in the profile is real work.
2. **The baselines improve where the stack used to bind and are flat elsewhere**: grvingt
   50 + 12 goes 35.7 -> **31.7 ns/point** against session 3 (+12 %), grdix 150 + 105 goes
   47.4-48.2 -> **46.0-47.3** and 127 + 128 stays at 39.7 (session 4: 37.1).
3. **The memory wall is a factor 2.4 away on grvingt and 5 away on grdix** (§3).  Each point
   is written to DRAM once and read once, which costs **48.7 measured bytes of DRAM traffic
   per 16-byte point on grvingt** and **~28 on grdix**; the router's baseline uses 41 % of
   grvingt's sustainable bandwidth and 19 % of grdix's.
4. **`192 + 63 --lossy` still collapses** and no longer recovers: 51 % dropped in every round
   (session 4 §2 saw it latch per round).  Lossless at the same split is clean.

## 1. Baselines, np 1

grvingt-11, 50 senders + 12 receivers (the reference split of sessions 1 and 3):

| mode | routed | ns/point | blocks/s | dropped |
|---|---|---|---|---|
| lossless | 1.6 G | 31.5 / 32.2 / 32.2 | 380-389 K | 0.000 % |
| lossy | 1.5-1.6 G | 31.4 / 31.8 / 31.8 | 385-390 K | 3.2-3.3 % |

grdix-5, three splits of the 255-worker team, lossless and lossy:

| split | mode | routed | ns/point | blocks/s | dropped |
|---|---|---|---|---|---|
| 150 + 105 | lossless | 3.2 G | 46.1 / 46.8 / 47.3 | 772-793 K | 0.000 % |
| 150 + 105 | lossy | 3.2-3.3 G | 45.8-46.1 | 794-799 K | 0.008-0.066 % |
| 192 + 63 | lossless | 2.5-2.8 G | 73.4 / 76.4 / 76.9 | 611-676 K | 0.000 % |
| 192 + 63 | lossy | **1.8 G** (push 3.7 G) | 52.1-52.4 | 890-895 K | **51 %** |
| 127 + 128 | lossless | 3.2 G | 39.7 / 39.7 / 39.8 | 777-779 K | 0.000 % |
| 127 + 128 | lossy | 3.2 G | 39.6-39.7 | 778-780 K | 0.005 % |

The senders' floor, for scale: 13.0 ns/point on grvingt, 10.4 on grdix (the `raw` line, PRNG
and the benchmark's destination pick, no push).  **Do not compare ns/point across the two
nodes**, only ratios.

## 2. Where the cycles go now

`perf record -F 499 --all-user` on one 3 s round, source lines summed by region of
`include/router/workers.hpp` (line ranges from the function boundaries), 99.8 % accounted:

| region | grvingt 50 + 12 | grdix 150 + 105 |
|---|---|---|
| **the free ring** (`free_push` / `free_pop_many` / `refill`) | **0.13 %** | **0.38 %** |
| `seal` + `stage_line` (the install wait) | 0.44 % | 1.23 % |
| `Router_Push` body | 23.00 % | 12.70 % |
| `Router_Grab` / `Router_Release` / `Router_Pop` | 3.80 % | 29.30 % |
| `memmove` (private line -> block, block -> out) | 5.42 % | 25.69 % |
| `tools.hpp` (PRNG, `murmur128`) | 34.82 % | 24.08 % |
| the benchmark's own lines | 26.82 % | 4.21 % |

Two things to read off it.  The free ring and the install wait together are **under 2 %** on
both machines, against 36 % (grdix, session 4 §4) and 41 % (grvingt np 2, session 1 §1)
before: the MPMC ring did what sharding was supposed to do.  And the two machines are now
limited by different halves of the router -- grvingt by the sender's push (23 %) and grdix by
the receiver's read-out plus the copies (55 % between `Grab`/`Pop` and `memmove`).

`router_bench.cpp:49`, the `x % F` destination pick, is **21.8 % of all samples on grvingt**
and 1.5 % on grdix.  Sessions 1 §3 and 3 §6 said exactly this: the divider is the single
biggest line on Intel and free on Zen.  It is still there, and it is still harness cost.

## 3. The memory wall

Every point is written into a block that is far larger than any cache by the time it is read,
so each point costs one DRAM write and one DRAM read, plus the read-for-ownership the write
pulls in.  That makes aggregate memory bandwidth a hard ceiling on the router, and it is worth
knowing how close it is.

### The machines

`membw`, an OpenMP STREAM in the session's scratch (`~/membw.c` on the Nancy home): arrays
first-touched in parallel, `OMP_PROC_BIND=spread OMP_PLACES=cores`, best of 5.  "DRAM" counts
the read-for-ownership that every write of a line not in cache costs, which is what the memory
controllers see; "STREAM" is the classic convention that counts only the bytes the program
moved.

| node | threads | read | copy | scale | triad |
|---|---|---|---|---|---|
| grvingt, 3 x 4 GiB | 32 cores | 161 | 190 | 190 | 189 |
| grvingt | 64 PU | 204 | 185 | 184 | 181 |
| grvingt, one socket | 16 cores | 81 | 97 | 97 | 96 |
| grdix, 3 x 8 GiB | 256 cores | 543 | 491 | 482 | 490 |
| grdix | 512 PU | 531 | 485 | 476 | 486 |
| grdix, one socket | 128 cores | 271 | 245 | 241 | 245 |

(GB/s of DRAM traffic.  Sockets scale exactly 2x on both.)

### Bytes of DRAM per point, measured

Not modelled: counted with the memory controllers' own counters, as the difference between a
1-round and a 3-round run so that setup and the `raw` phase cancel.

**grvingt**, `uncore_imc/cas_count_read/` + `cas_count_write/`, 50 + 12: 425352 MiB over
9.16 G points = **48.7 bytes per point** (27.3 written, 21.4 read).  The model -- 16 B of point
written, 16 B of read-for-ownership, 16 B read back -- predicts 48.  It is a coincidence worth
noting that the split is not the model's 2:1 read:write: some block reads are still being
served by the 44 MB of L3 (the inboxes hold ~50 MB in flight), and the writes carry more than
the points.

**grdix** exposes no named DRAM event; `amd_umc_*` has an empty `events/` directory, so the
data fabric's `local_processor_{read,write}_data_beats_cs0` was used and **calibrated against
`membw`**, whose traffic is known: 602 bytes of machine-wide traffic per channel-0 beat.  At
150 + 105 the router draws 0.0473 beats per point, i.e. **~28.5 bytes per point, +/- 25 %** --
consistent with 32 (write + read, no read-for-ownership: Zen combines full-line writes) and
clearly under grvingt's 48.7.

### How far the wall is

| | sustainable DRAM | bytes/point | ceiling | baseline | uses |
|---|---|---|---|---|---|
| grvingt, 50 + 12 | 190 GB/s | 48.7 | **3.9 G pts/s** | 1.58 G | 77 GB/s, **41 %** |
| grdix, 150 + 105 | 490 GB/s | ~28.5 | **~17 G pts/s** | 3.26 G | 93 GB/s, **19 %** |

So memory bandwidth does not bind either machine today, but on grvingt it is only **2.4x**
away and would be the next wall after the sender's push path.  On grdix there is **5x** of
headroom, which says the AMD node is limited by the router and not by its memory system.

A second, cruder reading of the same wall, for the sceptical: `membw`'s last block writes
16-byte points into a buffer and reads them back, which is the router's pattern with none of
the router in it.  It reaches **2.38 G pts/s on grvingt** and **8.6 G on grdix** -- 1.5x and
2.6x the baselines.  These are lower bounds on the machine, not the ceiling: the interleaved
16-byte format vectorizes badly, which is why they sit under the copy kernel's implied 3.9 G
and 17 G.

## 4. `192 + 63 --lossy` has stopped recovering

Session 4 §2 found this split bistable: round 0 nearly always collapsed and later rounds
climbed out to 2.6-2.7 G.  Now all three rounds collapse identically -- 1.8 G routed against
3.7 G pushed, **51 % dropped** (all `DROPPED_RECV`), service 78-83 % busy.  Lossless at the
same split is clean at 2.5-2.8 G with nothing dropped, and 150 + 105 lossy is clean too, so
the diagnosis stands: with drops returning blocks straight to the free ring the senders have
no back-pressure, and 3:1 is past the point where 63 receivers can drain them.  Relieving the
free stack has made this **worse**, exactly as session 4 §3 predicted it would: the stack was
the throttle.

## What to change

- **Give lossy mode real back-pressure.**  It is now the only thing standing between a 3:1
  split and throwing away half the points, and the free ring no longer throttles anything.
- **Fix the benchmark's `x % F`** (21.8 % of samples on grvingt).  It has been the top line on
  Intel in three sessions; a multiply-shift removes it.
- **grdix's receiver path is where its cycles are** (55 % between `Router_Grab`/`Pop` and
  `memmove`), and the bench pops one point at a time.  Measuring `Router_Grab` in place
  against `Router_Pop` would separate the router's cost from the benchmark's.
- Nothing to do about memory bandwidth yet, but note the number: **48.7 bytes of DRAM per
  16-byte point on Intel**.  Anything that makes a point cross DRAM twice would put grvingt
  straight into the wall.

## Reproducing

```bash
# baselines
mpirun -np 1 --bind-to none build/examples/router_bench --senders 50 --receivers 12          # grvingt
mpirun -np 1 --bind-to none build/examples/router_bench --senders 150 --receivers 105        # grdix
#   ... --lossy for the second half of the table

# the machine's bandwidth (membw.c is ~140 lines of OpenMP, in the session scratch)
gcc -O3 -march=native -fopenmp -o membw membw.c
OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=256 ./membw 8 5     # 3 arrays of 8 GiB

# bytes of DRAM per point: difference two runs so setup and the raw phase cancel
for r in 1 3; do perf stat -a -e uncore_imc/cas_count_read/,uncore_imc/cas_count_write/ -- \
    mpirun -np 1 --bind-to none build/examples/router_bench --senders 50 --receivers 12 --rounds $r; done
# on AMD there is no named DRAM event: use one channel and calibrate it against membw
EV=amd_df/local_processor_read_data_beats_cs0/,amd_df/local_processor_write_data_beats_cs0/

# the regions of the profile
perf record -F 499 --all-user -o g.data -- mpirun ... --rounds 1
perf report -i g.data --stdio --no-children -g none --sort srcline    # then bucket by workers.hpp line
```

---

# Session 6 -- grdix, where the cycles go at 150 + 105, and the 2.1x that was hiding there

Same day and same node as session 5 (**grdix-5**, 2 x EPYC 9754, 256 cores / 512 PU, 2 NUMA
nodes, 32 L3 domains of 16 MB), stock tree at `739b7ea` plus two out-of-tree patches, named
where used.  The first is a **role probe**: twelve lines in `RouterPlacement::report_measured`
printing `thread_cpu[]` by role, so that `perf -C` is an exact role filter (`perf record` also
needs `-a`, or the samples carry no CPU and `-C` silently returns nothing).

## Summary

| | stock | with non-temporal stores |
|---|---|---|
| 150 + 105 | 3.3 G pts/s, 46.0 ns | **6.9 G pts/s, 21.8 ns** |
| 192 + 63 | 2.5-2.8 G, 76.4 ns | **7.2 G, 21.2 ns** |
| 96 + 31, socket 0 only (128 cores) | 4.4 G, 22.0 ns | **6.2 G, 15.4 ns** |
| the sender's 8 KB copy | 33.8 % of sender cycles | 5.5 % |
| the sender's DRAM fills | 2.27 G (49 % from the far socket) | 0.04 G |
| the receiver's remote-cache fills | 1.53 G | 0.02 G |
| receiver IPC | 0.39 | 0.83 |

Three findings, in order of what they cost:

1. **The points were travelling cache-to-cache across the sockets, and both roles were paying
   for it** (§2, §3).  The receivers read almost nothing from DRAM: 1.53 G of their fills come
   from *the other socket's cache*, where the sender that wrote the block still holds the
   lines.  The senders, meanwhile, fetch every line they are about to overwrite (write
   allocation), half of those from far DRAM.
2. **One 8 KB `memcpy` per 512 points is a third of every sender's cycles** (§2), and it is
   slow because of that write allocation, not because of the copy.  Replacing it with
   non-temporal stores is worth **2.1x** on the whole node and takes the copy out of the
   profile (§4).  `ctest` passes all three.
3. **Cross-socket traffic is what stops the node scaling.**  The same 128-core team is
   4.4 G pts/s on one socket and 2.9 G split over two (§3) -- and with non-temporal stores one
   socket does 6.2 G while both do 7.2 G, so the second socket adds 16 % (§5).  This
   contradicts session 2 §3, which found socket crossing irrelevant; that was measured when
   the free-stack CAS was 74 % of sender cycles and masked everything.

## 1. The three roles, separately

`perf record -a -F 999 --all-user`, one 3 s round at 150 + 105 (3.3 G pts/s), samples filtered
by the role's CPUs, source lines of `include/router/workers.hpp`:

| senders (150 threads) | share of sender cycles |
|---|---|
| `__memmove_avx512_unaligned_erms` (the copy in `stage_line`) | **33.8 %** |
| `Router_Push` body, the two stores into the private line (`:164`, `:167`) | 15.7 % |
| the PRNG (`tools.hpp:130-143`) | 28.8 % |
| the benchmark's `x % F` (`router_bench.cpp:49`) | 2.3 % |

| receivers (105 threads) | share of receiver cycles |
|---|---|
| **`Router_Pop`'s two loads of the point out of the block** (`:245`, `:239`) | **84.0 %** |
| the benchmark's `murmur128` | 4.2 % |
| the inbox ring (`ring.hpp`), `free_push` (`:21`), `Router_Test_drained` (`:256`) | 2.7 % |

The service thread is one thread of 256 and idle (`ring.hpp:50` and `service.hpp:281`, its
sweep, are its own top lines; 38 % of its cycles are `MPI_Testsome` over an idle network, as
session 4 §3 already noted).

On the two receiver lines: 244 and 245 are the `.key` and `.val` halves of one 16-byte point,
i.e. **one cache line for four points**, and 239 is the loop-top test of the next iteration.
Cycle sampling on this machine is not precise (`cycles:pp` is unsupported, `-e ibs_op` would be
needed), so a load's cost lands on the instruction after it: what the three lines say together
is that essentially all of a receiver's time is the point load, not that the compare is
expensive.

## 2. What the `memmove` is, and why it is slow

It is one call per staged line, at `workers.hpp:132` in `stage_line`:

```c
memcpy(pts + (size_t) k * swc_linesize, line, (size_t) n * sizeof(Point));
```

`swc_linesize` is 512 points, so **8 KB from the sender's private write-combining line into
its slot of the destination block**.  glibc resolves it to `__memmove_avx512_unaligned_erms`
and, at 8 KB, takes the **4 x 64 B unrolled vector loop, not the ERMS `rep movsb` branch**:
85 % of the memmove's samples sit on four `vmovdqa64` stores of the loop body.

Those are ordinary stores, so each one that misses **fetches the line first** (write
allocation), and the profile is a stall on that fetch, not the copy.  Measured on the sender
CPUs over one round (`ls_any_fills_from_sys.*`, 9.8 G points routed, i.e. 2.45 G lines
written):

| fill source | count | share |
|---|---|---|
| near DRAM (own socket) | 1.156 G | 42 % |
| **far DRAM (other socket)** | **1.112 G** | **41 %** |
| remote cache | 0.357 G | 13 % |
| own L3 domain | 0.104 G | 4 % |

**2.73 G fills for 2.45 G lines written**: every line the sender writes is fetched from
somewhere first, and half of them come from the other socket.  Sender IPC is 0.97.

## 3. The receivers are not reading DRAM at all

Same measurement on the receiver CPUs, same round:

| fill source | count |
|---|---|
| **remote cache (the other socket's)** | **1.531 G** |
| own L3 domain | 0.055 G |
| near DRAM | 0.014 G |
| far DRAM | 0.005 G |

Of the 2.45 G lines the receivers read, 1.53 G arrive **from a cache on the other socket** and
19 M from DRAM.  The point stream is not going through memory: a sender writes a block into
its own cache hierarchy, and the receiver -- picked uniformly, so on the far socket half the
time -- probes it out across the fabric.  Receiver IPC is **0.39**, and 84 % of the cycles are
those two loads.

That is why session 5's DRAM figure was low (~28 bytes per point measured against a 48-byte
model): most of the traffic never reached DRAM.  It also predicts the socket test, at a
constant 128 cores and the same team shape:

| 96 senders + 31 receivers, 128 cores | routed | ns/point |
|---|---|---|
| `taskset -c 0-127`, one socket | **4.4 G** | **22.0** |
| `taskset -c 0-63,128-191`, 64 cores per socket | 2.9 G | 33.5 |

**1.5x for crossing sockets**, with the same thread count and the same code.

## 4. Non-temporal stores: 2.1x

The patch, in `stage_line`, in place of the `memcpy` (out-of-tree, AVX-512 only as written;
a shipping version wants a 256-bit path or a `memcpy` fallback):

```c
__m512i *dq = (__m512i *) (pts + (size_t) k * swc_linesize);
const __m512i *sq = (const __m512i *) line;
size_t nq = (size_t) n * sizeof(Point) / 64;
for (size_t q = 0; q < nq; q++)
        _mm512_stream_si512(dq + q, _mm512_load_si512(sq + q));
_mm_sfence();                  /* NT stores are not ordered by the release store below */
```

Both sides are 64-byte aligned by construction (`ROUTER_HDR_BYTES` is 64, the pool is
`aligned_alloc(64)`, a slot is `swc_linesize * 16` bytes and `swc_linesize` is a power of two
at least 4), so the aligned NT store is safe.  `ctest` passes `router_np1/2/4`.

| split | stock | NT stores |
|---|---|---|
| 192 + 63 | 2.5-2.8 G | **7.2 G** |
| 150 + 105 | 3.3 G | **6.9 G** |
| 127 + 128 | 3.2 G | 6.4 G |
| 105 + 150 | 3.0 G | 5.5 G |
| 85 + 170 | 3.2 G | 4.5 G |

Note the sign flip: stock wanted receivers (150 + 105 beat 192 + 63 by 30 %), and with the
fills gone the best split is the one with the **most senders**, because the receiver's read is
now cheap.  The fills say why: senders 2.73 G -> 0.24 G (DRAM 2.27 G -> 0.04 G), receivers'
remote-cache 1.53 G -> 0.02 G, replaced by 0.60 G of DRAM reads.  Receiver IPC goes 0.39 ->
0.83 for the same cycles, and the profile changes shape -- the copy falls from 33.8 % to 5.5 %
of sender cycles, and `Router_Push`'s two stores become the top sender lines.

**62 % of the receivers' new DRAM reads are far** (0.370 G against 0.226 G near), which is the
next thing to fix: the pool is first-touched in parallel by the whole team, so a block's pages
are wherever they landed.  Partitioning the pool by NUMA node and handing a sender a block on
its destination's node would make the receiver's read local.

## 5. And then the memory wall

With the fills gone, the one-socket configuration runs into DRAM, which is the wall session 5
was looking for.  Measured with the data fabric's channel-0 beats, calibrated against `membw`
on the same socket (511 bytes of machine traffic per beat, socket 0 alone):

| 96 + 31, socket 0, NT stores | value |
|---|---|
| routed | 6.22 G pts/s |
| DRAM per point | 34.7 B (the model says 32: one write, one read) |
| DRAM traffic | **216 GB/s** |
| that socket's copy bandwidth (`membw`) | 245 GB/s |
| | **88 % of the wall** |

The full node is not there yet: 7.2 G pts/s at 192 + 63 against 490 GB/s of bandwidth, i.e.
roughly half.  What it is short of is not bandwidth but the fabric -- one socket does 6.2 G
and two do 7.2 G.

## 6. On Intel the same patch is worth 3 %

grvingt-11, 50 + 12: 31.6 -> **30.9 ns/point**, 1.6 G either way.  The RFOs do disappear
(`offcore_requests.demand_rfo` 856 M -> 48 M, 18x), but the DRAM traffic per point does not
move (48.7 -> 52.9 bytes), so on that machine the fetches were being served by L3 rather than
by memory, and the copy was only 5.4 % of sender cycles to begin with.  Two open ends there:
its per-point traffic is ~1.6x the 32-byte model in both builds and I could not attribute the
excess, and glibc picks `__memmove_evex_unaligned_erms` rather than the AVX-512 variant.

## 7. Ruled out

- **Cache residency of the active destination blocks.**  105 destinations x 64 KB is 6.7 MB,
  which would fit an L3 domain, and shrinking it makes things **worse**, monotonically:
  `--block` 512 / 1024 / 2048 / 4096 / 16384 gives 2.1 / 2.2 / 2.3 / 3.3 / 3.5 G pts/s.  Small
  blocks cost more per block than they save in cache.
- **`--dests` below the receiver count** is refused outright (`dests_per_node must be at least
  the receivers per node`), so the fan-out cannot be narrowed to test residency that way.
- **The free ring**, again: 0.38 % of samples here, as in session 5.

## What to change

- **Stage lines with non-temporal stores** (§4).  2.1x on this machine, 3 % on Intel, nothing
  measured against it.  It needs an AVX2 path and a `memcpy` fallback, and the `sfence` before
  the release store is not optional.
- **Then partition the pool per NUMA node** and hand a sender blocks on its destination's node
  (§4): 62 % of the receivers' reads are far today.
- **Re-tune the default split after both.**  With the fills gone the optimum moves from
  150 + 105 to 192 + 63, i.e. towards senders.
- Session 4 §3's advice to shard the free ring is **done and no longer needed** (the MPMC ring
  did it), and its lossy back-pressure item is still open (session 5 §4).

## Reproducing

```bash
# the role probe: print thread_cpu[] by role at the end of RouterPlacement::report_measured,
# then filter samples by role -- and record system-wide, or -C returns nothing
perf record -a -F 999 --all-user -o a.data -- mpirun -np 1 --bind-to none build/examples/router_bench \
    --senders 150 --receivers 105 --rounds 1
perf report -i a.data --stdio --no-children -g none --sort srcline -C "$(cat cpu.senders)"

# where a role's cache lines come from
perf stat -a -C "$(cat cpu.senders)" -e cycles,instructions,ls_any_fills_from_sys.dram_io_near,\
ls_any_fills_from_sys.dram_io_far,ls_any_fills_from_sys.local_ccx,ls_any_fills_from_sys.remote_cache -- mpirun ...

# which memmove branch, at instruction level
perf report -i a.data --stdio --sort symbol -C "$(cat cpu.senders)" | grep memmove
perf annotate -i a.data --stdio --symbol=__memmove_avx512_unaligned_erms -C "$(cat cpu.senders)"

# crossing sockets, at a constant 128 cores
taskset -c 0-127        mpirun -np 1 --bind-to none build/examples/router_bench --senders 96 --receivers 31
taskset -c 0-63,128-191 mpirun -np 1 --bind-to none build/examples/router_bench --senders 96 --receivers 31

# Intel's RFOs
perf stat -a -e l2_rqsts.all_rfo,offcore_requests.demand_rfo -- mpirun ...
```

---

# Session 7 -- the non-temporal staging copy, shipped and A/B'd on both nodes

2026-09-08, **grdix-2** and **grvingt-12** (Grid'5000 Nancy), release builds of the committed
tree at `ba8b1eb` ("router: stage lines with non-temporal stores") against its parent
`6b26941` in a second worktree, same node, same flags, median of 3-4 rounds of 3 s.  Session 6
measured this change as an out-of-tree patch; this is the version in the tree, which differs
in shape only: `router_stream_copy` in `common.hpp`, with an AVX-512 path, an AVX2 path for
Rome-class machines and a `memcpy` fallback, and the `sfence` that orders the streaming stores
before the release store publishing the line.  `ctest` passes on both nodes; the engine's three
smoke tests pass.

## Summary

**It ships unconditional, because it wins in every configuration the engine can actually
run.**  The one regime where it loses is receiver-heavy, which the engine cannot reach: it
runs one dict thread per cache domain and puts everything else on producers.

| grdix-2, 256 cores, 32 L3 domains | before | after | |
|---|---|---|---|
| 220 + 32 (the engine's shape here) | 2.4 G pts/s, 191.4 ns | **6.9 G, 31.6 ns** | **2.9x** |
| 192 + 63 | 2.4 G, 79.8 ns | **7.1 G, 27.0 ns** | **3.0x** |
| 240 + 15 | 2.0 G, 189.6 ns | 5.2 G, 160.0 ns | 2.6x |
| 150 + 105 | 2.8 G, 54.0 ns | 6.6 G, 22.5 ns | 2.4x |
| 127 + 128 | 2.8 G, 45.5 ns | 6.2 G, 20.4 ns | 2.2x |

| grvingt-12, 32 cores, 2 L3 domains, one thread per core | before | after | |
|---|---|---|---|
| 29 + 2 (the engine's shape here) | 518 M pts/s, 55.9 ns | **558 M, 52.3 ns** | **+7.6 %** |
| 28 + 3 | 773 M, 36.4 ns | **837 M, 33.3 ns** | **+8.3 %** |
| 26 + 5 | 1.3 G, 20.6 ns | **1.4 G, 18.8 ns** | **+8.7 %** |
| 24 + 7 | 1.5 G, 16.0 ns | 1.4 G, 16.6 ns | **-3.8 %** |
| 20 + 11 | 1.2 G, 16.35 ns | 1.2 G, 17.2 ns | -5.0 % |
| 16 + 15 | 16.1 ns | 17.3 ns | **-7.2 %** |

Both machines are stable to under 1 % round to round at these sizes, so the small numbers are
real, not noise.

## 1. The crossover, and why it is only on the small node

On grvingt the sign flips at about **four senders per receiver**.  Below that the streaming
store is a loss of 4-7 %, above it a gain of 8-9 %.  The mechanism follows from what the two
machines do with a staged block:

- grvingt has **two L3 domains of 22 MB**, one per socket, and a destination's active block is
  64 KB.  With few senders the handoff often stays inside a socket's L3, so the write
  allocation the old `memcpy` paid was an L3 hit and the receiver's read was another.  Pushing
  the line to DRAM instead trades two cache hits for two DRAM accesses.
- grdix has **32 domains of 16 MB** over two sockets, and the destination is picked uniformly,
  so the reader is in the writer's domain about one time in 32.  Session 6 measured what that
  costs: 2.27 G DRAM fills per round on the sender side, 49 % of them from the far socket, and
  1.53 G cross-socket cache probes on the receiver side.  There is no cache hit to lose.

The engine's shape is what settles it.  `--dicts-per-node` is meant to be the number of cache
domains and never more than `--producers-per-node`, so a node runs **14 producers per dict
thread on grvingt and 7 on grdix** -- deep in the region where the streaming store wins on
both.  A receiver-heavy Router team is legal and pointless.

## 2. Lossy mode delivers nothing inside a node

Measured because the question was asked directly, and the answer is unambiguous.  `routed` is
the delivered rate; the printed ns/point is **per pushed point**, so a lossy run that throws a
quarter of its points away looks 20 % cheaper per point while delivering no more:

| grvingt 26 + 5 | routed | ns/point | dropped |
|---|---|---|---|
| lossless, before | 1.2 G | 21.0 | 0.000 % |
| lossy, before | 1.2 G | 16.5 | **23.8 %** |
| lossless, after | 1.4 G | 18.7 | 0.000 % |
| lossy, after | 1.4 G | 16.5 | **12.5 %** |

Same delivered rate, a quarter of the work wasted.  At 24 + 7 lossy is slightly *worse*
(1.4 G against 1.5 G lossless), and session 5 §1 has the extreme case: grdix 192 + 63 lossy
delivers 1.8 G with 51 % dropped where lossless is clean at 2.5-2.8 G.  Across sessions 3, 5
and 7 there is **not one single-node configuration where lossy delivers more than lossless**.

Why there is nothing to gain: the free ring is the only flow control the Router has.  In
lossless mode a sender that finds no block waits, which paces the senders to the receivers'
drain rate, keeps the queues short and the handoff cache-warm.  Dropping returns the block at
once, so senders never wait, the inboxes stay full, the write-to-read delay grows past cache
residency and the receivers slow down -- which drops more (session 4 §2).  Lossy earns its
place on the inter-node path, where the algorithm tolerates loss and the alternative is
blocking on the network.

**As a diagnostic it is useful but perturbing.**  The three drop counters do localise the
imbalance -- `svc` for no block installed, `net` for credit exhausted, `recv` for a full inbox
-- but what they report is the imbalance of a system with no flow control: 192 + 63 reads as
"51 % dropped at the receivers" while the lossless run at that split is clean.  The
non-perturbing form of the same signal is already counted and simply not printed:
`ROUTER_STALL_OUT` (a block parked because a receiver's inbox was full) and `ROUTER_STALL_IN`
(no free block to repost a receive).  Adding them to the round line, and a per-receiver spread
of delivered points beside the per-node push spread, would make a lossless run self-diagnosing;
the service thread already walks the per-thread counters to aggregate them, so only a small
API addition is missing.

## What to change

- Nothing about the copy: it is in the tree and verified on both machines.
- **Print the two stall counters** on the bench's round line, and stop reading single-node
  lossy runs as performance numbers (§2).
- Still open from session 6: **partition the pool per NUMA node** so a sender streams into
  memory local to its destination.  62 % of the receivers' DRAM reads were far after this
  change, and on one socket alone the router already runs at 88 % of that socket's bandwidth
  (session 6 §5), so page placement is what the next factor has to come from.

## Reproducing

```bash
# the A/B, on the node, against the parent commit -- the worktree registration lives in the
# shared .git, so give the path the node's name or `git worktree prune` first
git worktree add /tmp/pre-nt-$(hostname -s) 6b26941
cmake -S /tmp/pre-nt-$(hostname -s) -B /tmp/pre-nt-$(hostname -s)/build -DCMAKE_BUILD_TYPE=release
make -C /tmp/pre-nt-$(hostname -s)/build router_bench

for b in /tmp/pre-nt-$(hostname -s)/build build; do
    mpirun -np 1 --bind-to none $b/examples/router_bench --senders 220 --receivers 32 --rounds 3
done
```

---

# Session 8 -- gros, four hosts: the multi-host path runs at the wire, and lossless deadlocked

2026-09-08, **gros-68, gros-71, gros-90, gros-91** (Grid'5000 Nancy): one socket of Xeon Gold
5220 each, 18 cores / 36 PU, **one** 24.8 MB L3, 1 NUMA node, 92 GB, **25 GbE**.  Release build
of `4420082` (the streaming staging copy of session 7 included), worktree `~/mitm-gros` on
`bench-gros`.  One MPI rank per host, **14 senders + 3 receivers + 1 service = 18 threads on 18
cores**, no oversubscription.  Rounds of 3 s, the steady round of 2-3 reported.

## 0. Open MPI does not work between these hosts out of the box

Every node's global IPv6 address on `br0` is a **/128**, so the TCP BTL cannot conclude that two
nodes share a subnet: it builds an empty reachability graph and a rank aborts with
`MPI_ERR_INTERN` inside the first collective, printing "Unable to find reachable pairing between
local and remote interfaces".  IPv4 is a normal /20 and raw connectivity is fine (0.28 ms, TCP
connect ok).  What works:

```
--mca pml ^ucx --mca btl tcp,self --mca btl_tcp_if_include 172.16.64.0/20 \
--mca btl_tcp_disable_family 6 --prtemca oob_tcp_if_include br0 --prtemca plm_ssh_no_tree_spawn 1
```

Two traps.  `mpirun ... hostname` proves nothing -- `hostname` makes no MPI call, so it tests the
launcher and passes while every rank would abort on its first message.  And a filtered `mpirun`
hides the abort you need to read: two runs looked like silent failures because the pipeline
swallowed it.

## 1. One node, for reference

18 threads on 18 cores.  `raw` is the senders generating points without pushing, i.e. the
benchmark's own floor:

| senders + receivers | raw | routed | ns/point |
|---|---|---|---|
| 16 + 1 | 1.9 G pts/s | 300 M | 53.4 |
| 15 + 2 | 1.8 G | 598 M | 25.1 |
| **14 + 3** | **1.7 G** | **889 M** | **15.7** |
| 12 + 5 | 1.5 G | 841 M | 14.3 |
| 9 + 8 | 1.2 G | 614 M | 14.7 |

Raw is 120-133 M pts/s per sender thread (7.5-8.4 ns/point), so the Router's marginal cost at the
peak is 7.4 ns/point.  A receiver saturates at **~300 M pts/s**: one delivers 300 M, three deliver
889 M, and past three the senders are the constraint.

## 2. The wire, in the Router's own shape

`netbw`, one MPI thread per rank, non-blocking sends to every peer with a bounded window, 32
always-posted receives, 256 rotating send buffers so a message is read from cold memory as a block
is.  Three passes per size; the figure is the mean received rate per rank, with the spread across
ranks:

| message | received per rank | spread | aggregate |
|---|---|---|---|
| 8 KB | 0.58-0.62 GB/s | **0.00-0.82** | 2.3-2.5 GB/s |
| 32 KB | 0.72-0.99 | **0.00-1.19** | 2.9-4.0 |
| 65600 B (the Router's block) | 1.37-1.39 | 1.35-1.41 | 5.5 |
| 256 KB | 1.81-1.83 | 1.78-1.86 | 7.2-7.3 |
| 1 MB | 1.93-1.94 | 1.89-1.97 | 7.7-7.8 |

Monotone in the message size, saturating at **1.94 GB/s per host each way** from 1 MB up, which is
62 % of a 25 Gb/s link per direction, both directions at once.  In the Router's units that is
**121 M points/s per host** at 1 MB against 86 M at the default block size.

Two things to note about the small sizes.  The spread includes **0.00**: at 8 KB and 32 KB a rank
can be starved outright, so those rows are not just slow, they are unfair.  And the eager limit
here is 65536 bytes (`ompi_info --param btl tcp`), so the Router's 65600-byte block is 64 bytes
into the rendezvous path -- but rendezvous is *not* the problem: 65600 B beats every eager size
measured.

**A first version of this table was wrong and the mistake is worth recording.**  It counted a send
when it was *posted*, and an eager send completes into an internal buffer long before the data is
on the wire, so "sent" overstated by up to 50 % below the eager limit; and it printed rank 0's own
rate, which is the favoured rank.  Together they made 32 KB look faster than 65600 B, i.e. an
"eager cliff" that does not exist.  Count sends on completion, and report the spread.

## 3. Four hosts: the Router is at the wire

Uniform destinations, so 3 of 12 destinations are local and **three quarters of every node's
points cross the network** (the measurements confirm it: 116.3 M pushed, 87.2 M net).

| configuration | routed | pushed per host | net per host | egress | msgs/s |
|---|---|---|---|---|---|
| `--local-only` (no network) | **3.6 G pts/s** | 893 M | 0 | 0 | 0 |
| defaults, 65600 B messages | 465 M | 116 M | 87.2 M | 1.40 GB/s | 21.3 K |
| `--block 8192`, 128 KB | 550 M | 137 M | 103.1 M | 1.65 GB/s | 12.6 K |
| `--block 16384`, 256 KB | 598 M | 149 M | 112.1 M | 1.79 GB/s | 6.8 K |
| `--block 32768`, 512 KB | 610 M | 153 M | 114.4 M | 1.83 GB/s | 3.5 K |
| **`--block 65536`, 1 MB** | **625 M** | 156 M | 117.1 M | **1.87 GB/s** | 1.8 K |

Against §2's wire figures -- 1.38 GB/s at 65600 B and 1.94 at 1 MB -- the Router reaches **101 %
and 96 %**.  There is nothing to win in the code: the multi-host path is network-bound, and the
only lever is the message size, worth **+34 %** (465 -> 625 M pts/s).  The local-only control also
validates the whole setup: 893 M per host is the 889 M of §1.

**The network costs a factor of 5.8** on this cluster: 3.6 G local against 625 M once three
quarters of the points leave the node.  A 25 GbE cluster is a poor match for a router that moves
16 bytes per point at a billion points per second per node.

## 4. Lossy buys nothing here either

| configuration | routed | pushed per host | dropped |
|---|---|---|---|
| defaults, lossless | 465 M | 116 M | 0 |
| defaults, `--lossy` | 465 M | 449 M | **74 %** |
| `--block 65536`, lossless | 625 M | 156 M | 0 |
| `--lossy --credit 16 --block 65536` | 664 M | 426 M | **61 %** |

Same delivered rate at defaults for four times the pushing, and 6 % more in the tuned pair for
three points in five thrown away.  Session 7 §2 found the same on one node; the inter-node path
does not change the conclusion.  Lossy's value here is diagnostic: the drop counters name the
place that is full (`svc` no block installed, `net` credit or slots exhausted, `recv` a full
inbox).

## 5. The bug: lossless deadlocked with the credit raised

Every one of these hung in round 0, reproducibly, and the first one was left spinning for 16
minutes at 1800 % CPU before it was killed:

| configuration | before the fix |
|---|---|
| `--credit 16` | **hang** |
| `--credit 64` | **hang** |
| `--credit 16 --block 65536` | **hang** |
| `--credit 64 --block 65536` | **hang** |
| `--credit 16 --inbox 16` | **hang** |
| `--credit 4` (default), `--credit 8`, 10, 12, 14 | fine |
| any of the above with `--lossy` | fine |
| `--credit 16 --n-recv 64` | fine |
| `--credit 16` or 64 on **two** hosts | fine |

The state, from stacks of a hung rank: **all 14 senders inside `refill`** waiting on an empty free
ring, the three receivers spinning on **empty inboxes**, the service thread inside MPI progress
with nothing completing.  Every block in the pool was committed and none could come back.

The cycle: a node that has no block leaves its receive slots idle -- `repost()` takes a block from
the service's stash and, finding none, counts `ROUTER_STALL_IN` and returns.  A node that is not
receiving stalls every peer that sends to it, and a peer's send holds its block until it
completes, so those blocks never return either.  Raising the credit is what tips it: it lets more
of the pool sit in the network path at once, and `--inbox 16` tips it too because a full inbox
parks blocks instead of draining them.  Lossy never hangs because a dropped block goes straight
back through `release()`.  Two hosts never hang because one peer cannot commit enough.

**The census**, from a hung rank through gdb, is what settled it (14 senders, 3 receivers, the
default block, `--credit 16`):

| what | value |
|---|---|
| pool | 1277 blocks |
| the free ring (`free_in - free_out`) | **0** |
| the service's stash (`free_list`) | **0** |
| free send slots (`out_free`) | **0** of 32 |
| `credit_used` per peer | **{0, 16, 16, 0}** -- two peers at their full credit |
| parked, first four targets | none |
| failed reposts (`ROUTER_STALL_IN`) | **54,451,552** |
| blocks pushed into the three inboxes | 278, 278, 276 |

Every send slot in flight, both credit-limited peers full, the ring and the stash empty, and 54
million attempts to post a receive that found no block.  The pool had drained into the send path
of every node at once.

**The first fix was wrong and is worth recording.**  `c1a0e7d` gave the receive slots a reserve
that `release()` topped up before anything reached the senders.  It changed nothing: all four
configurations still hung.  The reason is in the census -- `release()` runs when a *send
completes*, which is exactly what is blocked, so the reserve is never fed.  Reverted in `97ab913`.

**The fix that works** (`6fb1c88`): a floor on the free ring.  `free_pop_many` takes a `floor`
argument; a sender leaves `n_recv` blocks behind, and the service passes zero, so the blocks a
sender may not have are exactly the ones its receives need.  What makes this one work where the
reserve did not is that a **local** receiver's `Router_Release` refills the ring without the
network being involved, so the floor is restored by traffic that cannot be blocked.  With one
node the floor is zero: nothing arrives, so a posted receive never needs replacing.

Measured after the fix, the four configurations that hung:

| configuration | routed |
|---|---|
| `--credit 16` | 480 M pts/s |
| `--credit 64` | 436 M |
| `--credit 16 --block 65536` | 501 M |
| `--credit 16 --inbox 16` | 463 M |

Note that raising the credit is still not worth anything once it works: 480 M against 459 M at
the default credit, and `--credit 16 --block 65536` is *worse* than `--block 65536` alone.  This
cluster's run-to-run spread is wide -- the same binary gave 429, 429, 460 and 459 M pts/s in four
passes of the same configuration, i.e. **+/- 7 %** -- so treat differences under 10 % here as
unresolved.

The floor's own cost, measured as a paired A/B against its parent commit, alternating the two
binaries three times over (defaults, the steady rounds of 3 s rounds):

| pass | with the floor | without |
|---|---|---|
| 1 | 461.4, 462.1 M pts/s | 459.9, 461.7 |
| 2 | 445.7, 445.9 | 453.8, 457.9 |
| 3 | 446.5, 447.5 | 458.5, 456.8 |

451 against 458 M on the means, i.e. **1.4 %**, in a harness where the same binary moves 3.5 %
between passes.  Call it free; it is certainly not worth a deadlock.

## What to change

- **Raise the default block for multi-node runs.**  4096 points is 65600 bytes and leaves a third
  of the network on the table; 65536 points (1 MB) is at the plateau on this fabric.  It is the
  only knob that matters across hosts.
- The engine will be network-bound on a cluster like this one, by a factor of ~6 against its
  local rate.  Anything that keeps points on the node -- a hash that favours local destinations,
  or a two-level routing scheme -- is worth more than any code tuning measured in these sessions.
- Sessions 2, 4 and 5 recommended `--credit 16` from measurements taken with **two ranks on one
  host**, where "the defaults cost 6.8x".  Across real hosts the same setting deadlocked until
  `c1a0e7d`, and it is worth only ~8 % once it works (465 -> 505 M at credit 10-12).  Shared-memory
  np 2 does not stand in for a network.

## Reproducing

```bash
# one rank per host, 18 threads on 18 cores, and the flags of §0 in $MCA
mpirun -np 4 --hostfile hostfile --map-by node --bind-to none $MCA \
    build/examples/router_bench --senders 14 --receivers 3 --rounds 3 --block 65536

# the wire, in the Router's shape (netbw.c is ~110 lines of MPI, in the session scratch)
mpirun -np 4 --hostfile hostfile --map-by node --bind-to none $MCA ./netbw 1048640 3 8 32 256

# the deadlock, before c1a0e7d: round 0 never ends, 18 threads spinning
mpirun -np 4 ... build/examples/router_bench --senders 14 --receivers 3 --credit 16
# where the threads are (the service is the one in MPI progress)
sudo-g5k gdb -p $(pgrep -x -n router_bench) -batch -ex "thread apply all bt 3"
```

---

# Session 9 -- 64 gros hosts: the Router reaches the fabric, once the receive slots scale

2026-09-08, **64 nodes of `gros`** (Grid'5000 Nancy), one MPI rank per host, 14 senders + 3
receivers + 1 service = 18 threads on 18 cores, no oversubscription.  Committed tree at
`6fb1c88` (the free-ring floor of session 8 included), release build in `~/mitm-gros`.  Uniform
destinations, so **63 of every 64 points cross the network**.  Rounds of 10 s, the second round
reported; `routed` is the aggregate delivered rate over all 64 hosts.

## Summary

**The default configuration runs at half the fabric; two knobs close the gap.**

| configuration | aggregate routed | per host egress | share of the wire |
|---|---|---|---|
| defaults (65600 B messages, 32 posted receives) | 4.0 G pts/s | 0.98 GB/s | at the wire for that size |
| 256 KB messages, 32 receives | 4.9 G | 1.21 GB/s | 73 % |
| 256 KB, 256 receives | 7.1 G | 1.74 GB/s | 105 % |
| 1 MB, 32 receives | 4.5 G | 1.11 GB/s | 58 % |
| 1 MB, 128 receives | 5.5 G | 1.36 GB/s | 71 % |
| **1 MB, 256 receives** | **7.9 G** | **1.94 GB/s** | **101 %** |
| 1 MB, 512 receives | 7.9 G | 1.94 GB/s | saturated |

**2x from the defaults**, and the tuned configuration sits exactly on the wire.  Nothing was
dropped in any lossless run.

## 1. The fabric first

`netbw`, one MPI thread per rank, 64 ranks, every rank sending to all others, window 8, 32
posted receives, 256 rotating send buffers:

| message | received per host | spread over 64 hosts | aggregate |
|---|---|---|---|
| 65600 B | 0.78 GB/s | 0.78-0.78 | 49.8 GB/s |
| 256 KB | 1.66 GB/s | 1.66-1.67 | 106.4 GB/s |
| 1 MB | 1.92 GB/s | 1.92-1.93 | 122.9 GB/s |

The fabric **scales**: per-host throughput at 64 hosts is within 1 % of the four-host figure at
1 MB (1.94), and the spread between the luckiest and unluckiest host is under 1 %, so no switch
stage saturates and nobody is starved.  What does not survive is the small message: 65600 bytes
falls from 1.38 GB/s at four hosts to **0.78** at 64.

## 2. Posted receives are the knob nobody had swept

At 256 KB, everything else default:

| posted receives | aggregate | per host |
|---|---|---|
| 8 | 3.6 G pts/s | 0.90 GB/s |
| **32 (the default)** | 4.9 G | 1.21 GB/s |
| 64 | 5.4 G | 1.32 GB/s |
| 128 | 5.9 G | 1.46 GB/s |
| 256 | 7.1 G | 1.74 GB/s |
| 512 | 7.1 G | 1.74 GB/s |

With 63 peers, 32 always-posted receives is too few: a peer's message cannot land unless a slot
is free, so its sender parks and stalls, and the effect compounds because a stalled sender's
node then sends less too.  The saturation point is **about four receives per peer** -- 256 for
63 peers, and at four hosts (3 peers) the default 32 was already past it, which is why sessions
5-8 saw nothing from this knob.  `n_recv` also sets the send slots, so both sides scale together.

The price is memory: the posted receives hold `n_recv * block_bytes`, i.e. 256 MB per node at
256 receives of 1 MB, and the pool is sized to cover them.

## 3. What does not help

- **Credit.**  `--credit 16` with 128 receives at 256 KB gives 4.9 G against 5.9 G at the
  default credit: raising it *costs* 17 %.  Sessions 2, 4 and 5 recommended 16 from
  shared-memory np 2 runs; on a real fabric it is worthless at best, and before `6fb1c88` it
  deadlocked.
- **Lossy.**  7.2 G delivered at 256 KB against 7.9 G for the best lossless run, while pushing
  615 M pts/s per host and **dropping 82 %** of it.  At 65600 bytes it is 5.0 G with 87 %
  dropped and the service thread at 60 % busy.  Four sessions on three machine classes now agree:
  lossy never delivers more.

## 4. The round's drain is not free

A 10 s push phase gives an 11.4 s round at 65600 bytes and 13.5 s at 256 KB, and the reported
rate divides by the whole round.  With 3 s rounds the same two configurations read 3.2 and 4.0 G
against 4.0 and 4.9 G at 10 s, i.e. **the drain costs a quarter of a short round**, and it grows
with the block size.  Measure multi-node throughput with rounds of at least 10 s.

## 5. What this cluster can and cannot do

| reference | aggregate | best measured | share |
|---|---|---|---|
| senders generating points, no routing (`raw`) | 106.7 G pts/s | 7.9 G | 7.4 % |
| one host's local routing x 64 (889 M each) | 56.9 G | 7.9 G | 13.9 % |
| the fabric at 1 MB | ~7.8 G | 7.9 G | at the wire |

A host manufactures 1.7 G points/s, which is 27 GB/s of point data, and its NIC carries 2 GB/s:
a **13-fold mismatch**.  With uniform destinations the machine therefore runs at a seventh of
its local capability, and no amount of tuning can change that -- the Router is already at the
wire.  What buys throughput on a fabric like this is **keeping points at home**: a sharding
scheme where most points resolve on the node that produced them is worth more than every
optimisation in sessions 1-9 put together.

## What to change

- **Scale `n_recv` with the peer count**, about four per peer, instead of the fixed 32; it is
  worth 2x at 64 nodes and nothing at 4, so a fixed default cannot serve both.  Cap it by the
  memory the posted receives cost.
- **Raise `block_points` for multi-node runs**: 65600-byte messages get 0.78 GB/s per host at
  this scale against 1.92 at 1 MB.
- Leave `credit` alone.
- Treat single-node lossy numbers, and now multi-node ones, as diagnostics only.

## Reproducing

```bash
# 64 hosts, one rank each; $MCA is the flag set of session 8 §0
mpirun -np 64 --hostfile hostfile.64 --map-by node --bind-to none $MCA \
    build/examples/router_bench --senders 14 --receivers 3 --seconds 10 --rounds 2 \
    --block 65536 --n-recv 256

# the fabric, same shape, for the ceiling
mpirun -np 64 --hostfile hostfile.64 --map-by node --bind-to none $MCA ./netbw 1048640 3 8 32 256
```

---

# Session 10 -- 32 grvingt hosts on Omnipath: 9.8 G points/s, and a transport that must be chosen

2026-09-08, **32 nodes of `grvingt`** (Grid'5000 Nancy), 2 x Xeon Gold 6130, 32 cores / 64 PU,
one MPI rank per host, **26 senders + 5 receivers + 1 service = 32 threads on 32 cores**.
Omnipath, port ACTIVE at 100 Gb/s.  `--seconds 20 --rounds 2`, the second round reported;
`routed` is the aggregate over all 32 hosts.  Tree at `534fb62`, rebuilt against the
`openmpi/4.1.6` module (`build416`).

## Summary

| | aggregate routed | per host egress | share of the wire |
|---|---|---|---|
| **256 KB messages, default 32 receives** | **9.8, 9.8, 9.9 G pts/s** | **4.76-4.78 GB/s** | **87 %** |
| 1 MB, default receives | 9.3, 9.1, 8.9 G | 4.33-4.50 GB/s | 81 % |
| 1 MB, 256 receives | 8.4, 9.2, 8.9 G | 4.06-4.45 GB/s | 78 % |
| 1 MB, 128 receives, **lossy** | 0.45 G | 0.09 GB/s | 99.1 % dropped |

Three repeats each, interleaved.  **9.8 G points/s from 32 hosts**, against 7.9 G from 64 gros
hosts on 25 GbE: this fabric is worth 4.6x per host (4.77 against 1.03 GB/s), and the Router
carries 87 % of it.

## 1. The transport has to be chosen, and the obvious choice is broken

`psm2` **does not work on this node image**, and it fails silently in a way that looks like
success.  The component loads and immediately unloads itself, so the `cm` PML has nothing to
open; with the TCP BTL still enabled the job then runs *over Ethernet* while the command line
says psm2.  Exclude TCP as well and the truth appears: "At least one pair of MPI processes are
unable to reach each other".  What the paths actually give, same 32 hosts, 256 KB messages:

| path | per host, all-to-all | `osu_bw` point to point |
|---|---|---|
| psm2 (`-mca mtl psm2 -mca pml ^ucx,ofi -mca btl ^ofi,openib`) | no connectivity | -- |
| TCP over 10 GbE (what the above silently becomes) | 1.04 GB/s | 1.39 GB/s |
| TCP over IPoIB (`ib0`) | 0.78 GB/s | -- |
| libfabric, `ofi` MTL | see §2 | **5.63 GB/s** |
| **UCX (`--mca pml ucx`)** | **5.49 GB/s** | 4.52 GB/s |

The two networks are physically separate with no L2 or L3 interconnection, so each figure
belongs to exactly one fabric: 1.04 GB/s is Ethernet, 0.78 is the Omnipath hardware driven
through its IP emulation, and 5.49 is the same hardware driven properly.  Two operational
notes: the system Open MPI (5.0.7) has no working Omnipath path at all, so
`module load openmpi/4.1.6` is required and the binaries must be rebuilt against it; and
`module load ... | tail` silently does nothing, because the pipe puts the module command in a
subshell and its PATH change is discarded.

## 2. `MPI_ANY_SOURCE` is a 500x portability trap

Over libfabric's `ofi` MTL -- the fastest transport by `osu_bw`, at 5.63 GB/s -- the Router
collapses to **14.5 M points/s aggregate, 28 messages/s per host, a 20 s round taking 119 s**.
The Router posts all of its receives with `MPI_ANY_SOURCE`, and that provider has no fast path
for it; UCX and the TCP BTL both do.  Nothing in the Router is wrong, and nothing about the
fabric is slow: the same code and the same hardware differ by **500x** between two MPI
transports.  Worth knowing before a production run picks a transport by accident.

## 3. The knobs do not matter here, and one earlier reading was noise

On gros over 25 GbE, posted receives were worth 2x and the message size 2.5x (session 9).  On
this fabric neither separates: 256 KB with the **default 32** receives is the best and the most
reproducible configuration measured, 1 MB is 7 % worse, and raising the receives to 256 or 512
does nothing.  An early single run at 128 receives read 7.6 G and a later one at 256 read 9.3 G,
which I first reported as a receive-slot effect; three interleaved repeats show it was
run-to-run spread.  The lesson is the measurement discipline, not the knob: single runs here
scatter by 15 %, repeats of the winning configuration by 0.5 %.

## 4. Two open items

- **The per-host imbalance is large and new.**  Hosts push between 37 and 407 M points/s in the
  same round, a 5 to 10-fold spread, against 1.4-fold on gros.  Since nothing is dropped and the
  aggregate is stable, the fast hosts are absorbing the slow hosts' share; what makes a host
  slow is unknown.
- **Receiver-heavy splits did not run**: 24 + 7 and 20 + 11 at 1 MB either timed out at 200 s or
  ended early, plausibly because 224 to 352 destinations at 1 MB make the pool 2-3 GB per node
  and startup exceeds the cap.  Needs a longer cap to be measured, not a fix.

## 5. Lossy is destructive on a fast fabric

1 MB blocks, 128 receives: senders push **1.5 G points/s per host**, the run delivers **445 M in
total** and **drops 99.1 %**, against 7.6 G delivered by the same configuration lossless.  With
160 destinations and 1 MB blocks most pushes find no installed block and the whole line goes.
Five sessions, four machine classes, two fabrics: lossy has never delivered more than lossless.

## Reproducing

```bash
module load openmpi/4.1.6                 # never pipe this command: the pipe loses its PATH edit
cd ~/mitm-grvingt && cmake -S . -B build416 -DCMAKE_BUILD_TYPE=release && make -C build416 router_bench
mpirun --hostfile ~/nodes.32 --map-by ppr:1:node --bind-to none --mca pml ucx \
    build416/examples/router_bench --senders 26 --receivers 5 --seconds 20 --rounds 2 --block 16384

# the wire, and the check that a transport is really being used
mpirun --hostfile ~/nodes.32 --map-by ppr:1:node --bind-to none --mca pml ucx ./netbw 262144 3 8 32 256
mpirun ... -mca mtl psm2 -mca pml ^ucx,ofi -mca btl ^ofi,openib,tcp ./netbw ...   # aborts: psm2 is dead here
```

---

# Session 11 -- the staging fence and the line publication, priced on grdix and grvingt

2026-09-09, **grdix-3** (2 x EPYC 9754, 256 cores / 512 PU, 2 NUMA nodes, 32 L3 of 16 MB) and
**grvingt-11** (2 x Xeon Gold 6130, 32 cores / 64 PU, 2 NUMA nodes, 2 L3 of 22 MB), one rank per
node, **np 1**, tree at `66d0648` built `-DCMAKE_BUILD_TYPE=release`, plus five out-of-tree `#ifdef`
patches that each switch one mechanism off.  `router_bench --seconds 3 --rounds 2`, three repeats
interleaved, median of the six rounds, reported as the node's push rate.  `--senders 150
--receivers 105` on grdix, `26 + 5` and `50 + 12` on grvingt.  `--swc` is pinned throughout, since
the fence rate is exactly what it sets.

## The question

Is the `_mm_sfence()` of `router_stream_copy` (`common.hpp:92`) worth removing -- by an epoch scheme
where a sender fences once in a while, publishes the value of a counter the service advances, and a
sealed block is complete once every sender has passed its tag -- and is `complete()`'s per-line scan
(`service.hpp:41`) worth replacing by one atomic counter per block?

## The variants

| | define | what is off | sound |
|---|---|---|---|
| a | -- | the shipped code | yes |
| b | `ROUTER_NO_SFENCE` | the fence | no |
| c | `ROUTER_NO_PUBLISH` | the fence, the `n_valid` store, `zero_valid`, `complete()`'s scan | no |
| d | `ROUTER_FENCE_EVERY_STORE` | *adds* a fence per 64 B store: the positive control | yes |
| e | `ROUTER_NO_COMPLETE_SCAN` | `complete()`'s scan alone | no |
| f | `ROUTER_NO_LINE_PUBLISH` | the `n_valid` store and the scan; the fence and `zero_valid` stay | no |

b, c, e and f ship blocks whose data may not have landed, so the points they deliver are garbage.
Nothing in the benchmark depends on the values and no counter it reports is affected, so the rates
are valid -- with one caveat that runs through everything below: a receiver may now read a block
whose streaming stores are still in flight, which can cost the *receiver* and therefore **biases
every unsound variant's gain downwards**.  The positive control d is sound.  Note also that c, e and
f all make `complete()` answer yes at once, so they remove not only work but the **wait** for a
sealed block's last line -- §4 is where that turns out to matter more than any instruction.

## Summary

| | grdix 150+105 | grvingt 26+5 | grvingt 50+12 |
|---|---|---|---|
| the rate at the auto line size | 6.8 G pts/s (`--swc 256`) | 1.30 G (512) | 1.64 G (512) |
| **the fence there (b)** | **+0.5 %** | **-0.2 %** | **-0.2 %** |
| everything (c) | +2.2 % | +6.3 % | -3.1 % |
| the fence at `--swc 64` (b) | +5.8 % | +17.2 % | +11.0 % |
| everything at `--swc 32` (c) | +20.0 % | +21.7 % | +16.3 % |
| one fence per 64 B store (d) | **-68 %** | **-72 %** | -- |

**The fence is free at the rate it is fired.**  An `sfence` after a run of non-temporal stores costs
~100 ns wherever its latency is exposed (§3), and the shipped configuration fires one per staged
line, i.e. once per 2-8 KB, once per 2.5-5.7 us of a sender's time.  Sixty-four times that rate
costs 68-72 % of the node (d), so the experiment resolves fences perfectly well; one per line is
simply 1/64th of a mechanism.  What the c column contains is not ordering: on grdix it is
`zero_valid`'s byte stores (§4), on grvingt it is the completion wait (§4), and on both the fence's
share of it is nil.

## 1. The fence, line size by line size

grdix 150+105, push M/s, median of 6 rounds, spread +-2 %:

| `--swc` | lines/block | a | b (no fence) | c (all off) | e (no scan) |
|---|---|---|---|---|---|
| 32 | 128 | 4838 | +1.1 % | +20.0 % | +0.8 % |
| 64 | 64 | 6357 | +5.8 % | +10.2 % | -0.6 % |
| 128 | 32 | 6796 | -1.4 % | +0.6 % | +0.2 % |
| **256 (auto)** | 16 | **6769** | **+0.5 %** | +2.2 % | +0.5 % |
| 512 | 8 | 6784 | +1.4 % | +3.6 % | +0.1 % |

grvingt 26+5 (32 threads on 32 cores) and 50+12 (63 threads on 64 PU), the same:

| `--swc` | 26+5 a | b | c | 50+12 a | b | c |
|---|---|---|---|---|---|---|
| 32 | 781 | +18.4 % | +21.7 % | 1163 | +15.0 % | +16.3 % |
| 64 | 1067 | +17.2 % | +27.1 % | 1401 | +11.0 % | +12.4 % |
| 128 | 1309 | -0.7 % | +6.4 % | 1530 | +6.2 % | +5.0 % |
| 256 | 1306 | -0.1 % | +7.4 % | 1577 | +3.3 % | -0.6 % |
| **512 (auto)** | **1300** | **-0.2 %** | +6.3 % | **1635** | **-0.2 %** | -3.1 % |

On grvingt 50+12 the fence's cost follows a clean `1/swc` law: 15.0, 11.0, 6.2, 3.3 % is 207, 252,
260, 269 ns per staged line -- **one number, ~250 ns per fence in situ**, more than the ~100 ns of a
bare streaming loop (§3) because the fence also drains the sender's `dest[]` `fetch_add` and its
`n_valid` store.  On grdix the same fence never costs more than 6 %, and nothing at all at the auto
line size: 150 streaming senders keep the memory system busy while any one of them waits (§3).

## 2. What an `sfence` costs, and when the cost is exposed

`~/nt_fence.c`, out of tree: T threads each streaming `line` bytes at a time with
`_mm512_stream_si512` into its own 64 MB region, one `_mm_sfence()` every k lines.  ns per line:

| | grdix, 4 KB | grdix, 512 B | grvingt, 4 KB | grvingt, 512 B |
|---|---|---|---|---|
| 1 thread, no fence | 143 | 17.9 | 573 | 70.1 |
| 1 thread, fence per line | 240 | 119.6 | 705 | 180.3 |
| all cores, no fence | 1348 (150 thr) | 168.5 | 872 (32 thr) | 109.0 |
| all cores, fence per line | 1368 | 170.4 | 940 | 197.7 |
| the machine's write ceiling | 455 GB/s | 455 GB/s | 150 GB/s | 150 GB/s |

A fence costs **~100 ns** wherever its latency is exposed, on both machines and at both line sizes:
it is one write's trip to the coherence point and does not depend on how much was written before
it.  Whether the *machine* loses that time depends on how many streams there are to cover it.
grdix at 150 threads pays **1 %** for a fence per 512 B line; the same machine at 64 threads pays
40 %; grvingt, with 32 cores, pays 45 % and cannot hide it at all.  That is why the Router's own
numbers differ between the two machines, and why the fence is free on the bigger one.  A fence every
8 lines is already free everywhere (30.1 ns against 17.9 at one thread, 34.2 against 23.8 at 16).

## 3. `complete()`'s scan is free, so the atomic counter buys nothing

Variant e -- `complete()` answers yes without reading a single valid byte -- is within noise at every
line size on grdix (+0.8, -0.6, +0.2, +0.5, +0.1 %), **including `--swc 32`, where it skips 128
bytes per block and the service handles 1.3 M blocks/s**.  The scan is not on anybody's critical
path; the service thread has had slack in every configuration measured since session 4.  On grvingt
e does move the rate (+5.9 % at `--swc 512`), but §4 shows that is the wait it also removes, not the
scan.

So the counter's stated benefit is worth zero, and its cost -- one contended RMW per line where
there is a plain byte store today -- is not.  It keeps one real benefit, which is that it makes
`zero_valid` one store instead of L (§4).

## 4. What the c column really contains, and it is different on the two machines

grdix, variant f (no `n_valid` store, no scan; the fence and `zero_valid` kept):

| `--swc` | a | f | for comparison: b | e | c |
|---|---|---|---|---|---|
| 32 | 4820 | **+6.0 %** | +1.1 % | +0.8 % | +20.0 % |
| 256 | 6830 | **+0.3 %** | +0.5 % | +0.5 % | +2.2 % |

At the auto line size nothing here is resolvable: every variant is inside the +-2 % spread.  At
`--swc 32` the per-line store is +5.2 % (f less e) and the components do not add up -- 1.1 + 5.2 +
0.8 against 20.0 for all of them together -- so the small-line regime is superadditive, and it is
the only one where a variant went bimodal (b's six rounds spread 2650..5051 M/s).  What is in c and
in none of b, e, f is `zero_valid`'s **L single-byte atomic stores per installed block**
(`workers.hpp:98`), a loop no compiler can vectorise, run for all 32 blocks of a sender's refill:
about 13 % at L = 128, and lost in the noise at L = 16.

grvingt tells a different story, and it is the more interesting one.  At `--swc 512` (L = 8, the
auto value, 5 destinations):

| variant | push M/s | vs a |
|---|---|---|
| a | 1300 | -- |
| b (no fence) | 1298 | -0.2 % |
| e (no scan) | 1375 | +5.9 % |
| f (no store, no scan) | 1391 | +6.5 % |
| c (all off) | 1381 | +6.3 % |

Every variant that makes `complete()` answer yes gains the same ~6 %, and the one that does not
gains nothing.  With L = 8 there are 8 byte stores per 4096 points, far too few to be worth 6 %, so
what those three variants have in common is the only remaining candidate: **they never wait for a
sealed block's last line.**  With 5 destinations, one block per destination installed and ~10 us for
a sender to fill a 512-point line, a destination's throughput is bounded by one block per laggard
line, and the pending list is the bottleneck.  grdix, with 105 destinations, shows none of it.

This is the finding that matters for the epoch scheme: **the wait, not the ordering, is what the
completion protocol costs**, and an epoch scheme replaces a wait for one block's own lines by a wait
for every sender on the node to reach a checkpoint.  It moves in the wrong direction.

## 5. The senders bind, and the private line is worth more than everything above

grvingt, `--swc 256`, the same 31 worker threads split differently:

| split | push M/s |
|---|---|
| 8 + 23 | 424 |
| 14 + 17 | 748 |
| 20 + 11 | 1092 |
| 26 + 5 | 1310 |

More senders is more throughput all the way to 26 + 5, so the ~1300 M/s plateau of §1 is a
sender-side limit and not the receivers -- and, per §4, it is the completion wait behind it.

The knob that dwarfs every effect in this session is the private line itself.  grdix, 150+105, with
`--dests 1024` to fake the fan-out of 64 hosts x 16 dict threads:

| `--swc` | push M/s |
|---|---|
| **32 (what the auto rule picks at F = 1024)** | **5340** |
| 64 | 6194 |
| **128** | **6258** |
| 256 (private lines 4 MB per sender) | 5296 |

The auto rule (`connect.hpp:201`) gives the private lines a 512 KB budget, which at F = 1024 forces
`--swc 32` and **leaves 17 % on the table**; the optimum is 128, a 2 MB budget.  Session 2's
"`--swc` is nearly irrelevant above ~128" was measured at F = 3.

## What to change

- **Nothing about the fence.**  One instruction per 2-8 KB of staged points, worth between -0.2 %
  and +0.5 % at the shipped line size on three configurations of two machines, and ~100 ns when
  fully exposed.  The epoch scheme would trade it for a global grace period, a rate limiter on the
  counter, a livelock to design against in `refill`, the loss of a locally checkable invariant --
  and a *longer* completion wait, which §4 says is the one thing in this mechanism that costs real
  throughput.
- **Nothing about `complete()`'s scan** (§3): free at every L measured, on a thread with slack.
- **`zero_valid` deserves a patch** (§4) and needs no protocol change: its L atomic byte stores can
  be one `memset` of the block's valid bytes -- the block is the sender's own, off the free ring,
  until the install's release store publishes it -- or one counter per block.  Worth ~13 % at
  L = 128 and nothing at L = 16, so it is insurance for the high-fan-out case, not a win today.
- **Attack the completion wait at small fan-out** (§4): with F destinations there is exactly one
  installed block each, so a laggard line stalls a whole destination.  More blocks in flight per
  destination (a second installed block, or `--dests` above the receiver count) is the direction;
  it is worth ~6 % at F = 5 and nothing at F = 105.
- **Raise the private-line budget** (§5): 2 MB instead of 512 KB in `connect.hpp:201` changes
  nothing below F = 128 and is worth 17 % at F = 1024.  Re-measure with `--block 262144` at the same
  fan-out before settling it, since session 9 wants big blocks on a fabric.
- Keep the `--swc` sweep in the tool box: it moves the fence rate, the reservation rate and the
  private-line footprint at once, and four of this session's five findings came out of it.

## Reproducing

```bash
# the A/B patches: each is a ~5-line #ifdef in common.hpp / workers.hpp / service.hpp
cmake -S . -B build-b -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS=-DROUTER_NO_SFENCE
make -C build-b router_bench -j
for rep in 1 2 3; do for swc in 32 64 128 256 512; do for v in a b c; do
    mpirun -np 1 --bind-to none build-$v/examples/router_bench --senders 150 --receivers 105 \
        --swc $swc --seconds 3 --rounds 2 --quiet
done; done; done

# the fence on its own, the split probe, the fan-out probe
cc -O3 -march=native -fopenmp -o nt_fence nt_fence.c && ./nt_fence 150 512 1 2
mpirun -np 1 --bind-to none build-a/examples/router_bench --senders 20 --receivers 11 --swc 256 ...
mpirun -np 1 --bind-to none build-a/examples/router_bench --senders 150 --receivers 105 \
    --dests 1024 --swc 128 --seconds 3 --rounds 2 --quiet
```

# Session 12 -- VT-d halves Omni-Path: debian11 turns the IOMMU on and nobody said so

2026-09-11, **grvingt-6** and **grvingt-7** (2 x Xeon Gold 6130, 32 cores / 64 PU, 2 NUMA nodes),
`oarsub -t deploy -p grvingt -l host=2`, kadeploy'ed **debian11**, kernel `5.10.0-38-amd64`
(Debian 5.10.249-1, 2026-02-10), Open MPI 4.1.0, libpsm2 11.2.185, Omni-Path 100 Gb/s.
Not the Router: plain `flood.c`, a 2-rank `MPI_Send`/`MPI_Recv` ping-pong of one 400 MB buffer,
`mpiexec -n 2 --map-by ppr:1:node`.

## The question

`flood.c` carries a 2019 reference line, `grvingt-12 ... débit : 12101.5 Mo/s`, which is 97% of the
100 Gb/s wire.  The same command on grvingt-6/7 now gives **6371 Mo/s**.  Nothing in the program
changed.  What did?

## It is not the fabric, and it is not the transport

| checked | reading |
| --- | --- |
| OPA link | `Active`, `LinkSpeed 25Gb`, `LinkWidth 4`, **Link Quality 5 (Excellent)**, both ends |
| OPA error counters | every one **zero** on both ends; only the data/packet counters move |
| PCIe of `hfi1` `0000:5e:00.0` | `LnkCap 8GT/s x16`, `LnkSta 8GT/s (ok) x16 (ok)` -- gen3 x16, no downgrade |
| HFI NUMA node | 0, on both nodes |
| governor | `performance`, `intel_pstate` active, `no_turbo=0` |
| THP | `[always]` |
| transport | `--mca pml cm --mca mtl psm2` gives **6345**, i.e. the same as the default: PSM2 *is* carrying it |
| the TCP fallback, for scale | `--mca pml ob1 --mca btl tcp,self,vader` -> **1149** Mo/s.  Never the case here. |
| binding | default already lands on the HFI's socket: `numactl --cpunodebind=0` 7800, `=1` 5600 Mo/s |

The default run prints `UCX ERROR ibv_create_cq(cqe=4096) failed: Operation not supported` twice.
That is cosmetic -- `hfi1` has no full verbs CQ, UCX bails, OMPI falls back to psm2 and gets the
identical number.  `--mca pml ^ucx` silences it.  `--mca pml ob1 --mca btl ofi` crashes in `libucs`.

## Two ranks are not the measurement

`flood.c` does exactly one round trip, so it prices the cold pinning of the buffer along with the
transfer.  An iterated version (`flood2`, same ping-pong N times on the same buffer) separates them:

| iteration | Mo/s |
| --- | --- |
| 0 (cold) | 6324 |
| 1..9 (steady) | ~7900 |

So there are two losses, not one: ~20% of cold setup on top of a steady state that is itself only
64% of the wire.  `/usr/bin/time` over 20 iterations: `3.81 real 1.16 user 1.42 sys` -- 1.4 s of
kernel time for 2.1 s of transfer.

## The cause

Debian's bullseye kernel is built

```
CONFIG_INTEL_IOMMU=y
CONFIG_INTEL_IOMMU_DEFAULT_ON_INTGPU_OFF=y
# CONFIG_IOMMU_DEFAULT_PASSTHROUGH is not set
```

so **VT-d DMA remapping is on by default with no `intel_iommu=` anywhere on the command line**.
`/proc/cmdline` is just `root=UUID=... console=... modprobe.blacklist=nouveau`, and yet
`DMAR: Intel(R) Virtualization Technology for Directed I/O`, four DMAR units on queued
invalidation, and the HFI in `/sys/kernel/iommu_groups/77` with `type = DMA` -- a translated
domain, not identity.  Every `hfi1` SDMA descriptor and every TID receive therefore walks IOMMU
page tables and pays IOTLB invalidations.  Intel's own tuning guide asks for this to be off:
"Intel(R) Omni-Path Fabric Performance Tuning User Guide", Doc. No. H93143 Rev. 19.0 (April 2020),
section 3.3 "Do Not Enable intel_iommu" ("Setting intel_iommu=on in the kernel command line can
hurt verbs and MPI performance"), and Table 3's "Intel(R) VT for Directed I/O (VT-d): Disabled",
which section 2.2 extends to Xeon Scalable.  Two caveats on that citation: section 3.3 was *added*
in Rev. 7.0 (April 2017), not present from the start, and Table 5 recommends VT-d **enabled** on
Xeon Phi x200, where X2APIC needs it -- so it is not a blanket rule.  And the guide only
anticipates an operator adding `intel_iommu=on` and says to remove it from grub; it does not
cover a kernel that turns it on by default with nothing on the command line.

The 2019 reference was taken as an ordinary user on the standard environment of the day
(Debian 9/10), whose kernel defaulted it off.  The fabric never changed; the environment did.

## Proof

`kexec` into the *same* kernel and initrd with `intel_iommu=off` appended, both nodes, nothing
else touched:

| | cold (`flood`) | steady (`flood2`) |
| --- | --- | --- |
| as deployed, VT-d translated | 6371 | ~7900 |
| `intel_iommu=off` | **11636** | **12343** |
| 2019 reference | 12101 | -- |

12343 Mo/s is 99% of the payload rate of a 100 Gb/s link.  Both losses go away together: the cold
penalty was the IOMMU mapping of the buffer, the steady one the per-descriptor translation.
So **1.83x on the number `flood.c` prints**, 1.56x in steady state.

The other suspect was the mitigation stack -- this kernel reports `retbleed: Mitigation: IBRS`
(expensive on Skylake), `meltdown: PTI`, and `Clear CPU buffers` for mds / tsx_async_abort /
mmio_stale_data, none of which existed in 2019.  The IOMMU alone accounts for the whole gap, so
they cost little here and `mitigations=off` was not needed.

## What to do

Boot the compute nodes with `intel_iommu=off`.  After kadeploy, **once** per node:

```bash
kexec -l /boot/vmlinuz-$(uname -r) --initrd=/boot/initrd.img-$(uname -r) \
      --append="$(cat /proc/cmdline) intel_iommu=off" && systemctl kexec
```

~35 s per node, it skips the Dell POST and the deployed filesystem is untouched.  For something
that survives, put `intel_iommu=off` in `boot: kernel_params:` of a custom environment
(`kaenv3 -p debian11-min > env.yaml`, edit, `kaenv3 -a env.yaml`) and deploy that.

`iommu=pt` boots and does the right thing on this kernel -- grvingt-6 came up with
`iommu: Default domain type: Passthrough (set via kernel command line)` and group 77 reading
`identity` -- but its **bandwidth is unmeasured**, because grvingt-7 did not come back from that
second kexec and the walltime ran out.  Do not chain kexecs: the second one, from an already
kexec'ed kernel with `hfi1` loaded, is what hung the node; a hard `kareboot3` clears it.

This is not a `flood.c` curiosity.  It is the transport under every multi-node Router run on
grvingt, and sessions 10 and 11 were measured on nodes that were paying it.  **This needs root,
i.e. a deploy job.**  On a standard-environment node, check `dmesg | grep 'Default domain type'`;
if it says `Translated` the only routes are a deploy job with a custom `kaenv3` environment, or
asking the Nancy admins to add the flag to the standard kernel parameters.

## Reproducing

```bash
oarsub -t deploy -p grvingt -l host=2,walltime=2 sleep infinity
kadeploy3 -f $OAR_NODEFILE -e debian11-min -k
mpicc -O2 -o flood flood.c   # on both nodes, same path

mpiexec -n 2 --hostfile $OAR_NODEFILE --map-by ppr:1:node --allow-run-as-root ./flood 400000000
# 6371 Mo/s

# both nodes, then wait ~35 s
kexec -l /boot/vmlinuz-$(uname -r) --initrd=/boot/initrd.img-$(uname -r) \
      --append="$(cat /proc/cmdline) intel_iommu=off" && systemctl kexec

mpiexec -n 2 --hostfile $OAR_NODEFILE --map-by ppr:1:node --allow-run-as-root ./flood 400000000
# 11636 Mo/s
```

# Session 13 -- debian13 on grvingt: the IOMMU default got cheaper, and PSM2 stopped working under it

2026-09-11, **grvingt-4** and **grvingt-5** (2 x Xeon Gold 6130), `oarsub -t deploy -p grvingt
-l host=2`, `kadeploy3 -e debian13-big`, kernel **6.12.107+deb13-amd64**.  MPI is not in the
image: `module load openmpi` gives a Guix-built Open MPI **4.1.6** under `/gnu/store`, linked
against psm2 12.0 (the system `libpsm2` is 11.2.185).  Same `flood` / `flood2` 400 MB ping-pong
as session 12, `mpiexec -n 2 --map-by ppr:1:node`.  Four boots, `kexec` into the same kernel with
one flag changed each time.

## Two things have to be fixed before anything runs

**1. The OPA device is not called `hfi1_0` any more.**  debian13's rdma-core gives it a
persistent name, `opap94s0`.  `libpsm2` finds its units by scanning
`/sys/class/infiniband/hfi1_<N>` (`hfi_get_num_units()`), so it counts zero, and
`--mca pml cm --mca mtl psm2` fails with `PML cm cannot be selected`; Open MPI silently falls
back to the openib BTL at a third of the wire.  Fix, after every boot, on **every** node:

```bash
ip link set ib0 down; rdma dev set opap94s0 name hfi1_0; ip link set ib0 up
```

**2. Rename before the first MPI run, not after.**  Renaming a device that the openib BTL has
already opened leaves it wedged: psm2 *and* openib then hang until the node is rebooted.  That
cost an hour here.

## The four variants

Same procedure each time: fresh `kexec`, rename, wait for `PortState: Active` on both ends, run.
`psm2` is `--mca pml cm --mca mtl psm2`, `openib` is `--mca pml ob1 --mca btl openib,self,vader`,
`ucx` is `--mca pml ucx`.  The openib and ucx columns were taken in a separate pass with the
device left at its `opap94s0` name, so psm2 was out of the picture there.

| flag | domain type | group | psm2 cold | psm2 steady | openib cold | openib steady | ucx cold |
| --- | --- | --- | --- | --- | --- | --- | --- |
| *(none -- the default)* | Translated | `DMA-FQ` | **hangs** | -- | 4506 4533 4536 | ~4630 | 2567 2651 |
| `iommu.strict=1` | Translated | `DMA` | **hangs** | -- | 1609 1618 1629 | ~1605 | 1018 1032 |
| `intel_iommu=off` | -- (DMAR 0) | -- | 11729 11705 11714 | ~10500 | 4986 5045 5066 | ~5390 | 4459 4563 |
| `iommu.passthrough=1` | Passthrough | `identity` | 11513 11666 11739 | ~10800 | 5045 5064 5144 | ~5375 | 4056 4556 |

## Three results

**The kernel default moved from strict to lazy, and that is worth 2.9x.**  debian13 still builds
`CONFIG_INTEL_IOMMU_DEFAULT_ON_INTGPU_OFF=y`, so VT-d is still on with nothing on the command
line -- but 6.12 also has `CONFIG_IOMMU_DEFAULT_DMA_LAZY=y`, and the HFI's group now reads
`DMA-FQ` (deferred invalidation behind a flush queue) where debian11's 5.10 read `DMA` (an
invalidation per unmap).  Forcing `iommu.strict=1` reproduces the old behaviour exactly, and it
costs 2.9x on openib (4630 -> 1605) and 4.5x on ucx.  **That, not VT-d as such, is what made
session 12's debian11 numbers so bad.**

**But PSM2 does not work at all while the IOMMU translates.**  Under both `DMA-FQ` and `DMA` the
psm2 MTL hangs: the link is `Active` with LIDs assigned at both ends, `psm2_ep_num_devunits()` now
finds the unit, and the transfer simply never completes -- killed at 60 s, and given 90 s x3 plus
150 s in an earlier pass, for a transfer that takes 0.07 s when it works.  Three independent boots,
same result.  With `intel_iommu=off` or `iommu.passthrough=1` it works first try.  On debian11 /
5.10 psm2 ran fine under translated-strict, just slowly (session 12), so this is new with 6.12.
The practical consequence is blunt: **on debian13, Omni-Path either runs with the IOMMU out of the
way or it does not run.**

**`iommu.passthrough=1` is as fast as `intel_iommu=off`,** to within the noise, on every transport
(psm2 11513-11739 vs 11705-11729 cold, openib ~5375 vs ~5390 steady).  It keeps interrupt
remapping, so it is the better of the two and closes the item left open in session 12.

## What to do

Boot the grvingt nodes `iommu.passthrough=1`, and rename the OPA device before the first MPI run.
After kadeploy, once per node:

```bash
kexec -l /boot/vmlinuz-$(uname -r) --initrd=/boot/initrd.img-$(uname -r) \
      --append="$(cat /proc/cmdline) iommu.passthrough=1" && systemctl kexec   # ~50 s
# then, on every node, before any MPI:
ip link set ib0 down; rdma dev set opap94s0 name hfi1_0; ip link set ib0 up
```

Durable version: `iommu.passthrough=1` in `boot: kernel_params` of a custom `kaenv3` environment.

## One open item

psm2 steady state is ~10500-10800 here against 12343 on debian11 in session 12, while the cold
single round trip is the same (11.7 G against 11.6 G).  Different Open MPI (Guix 4.1.6 + psm2 12.0
against Debian 4.1.0 + psm2 11.2.185) is the obvious suspect, but it was not chased.  Not chased
either: whether the psm2 hang is hfi1's DMA mapping or its TID cache, and whether a newer
`libpsm2` or `hfi1` fixes it.

## Reproducing

```bash
oarsub -t deploy -p grvingt -l host=2,walltime=2 sleep infinity
kadeploy3 -f $OAR_NODEFILE -e debian13-big -k
module load openmpi          # Guix Open MPI 4.1.6, absolute path under /gnu/store
mpicc -O2 -o flood flood.c   # on both nodes, same path

# on both nodes, per boot, BEFORE any MPI run
ip link set ib0 down; rdma dev set opap94s0 name hfi1_0; ip link set ib0 up

mpiexec -n 2 --hostfile $OAR_NODEFILE --map-by ppr:1:node --allow-run-as-root \
        --mca pml cm --mca mtl psm2 ./flood 400000000     # hangs on the default boot
# kexec with iommu.passthrough=1, redo the rename, run again -> 11.5-11.7 G
```

# Session 17 -- grvingt: the prefetch is worth +50 % per dict thread, FILL gets +33 % from one too, and the team's shape decides the rest

2026-09-14, **grvingt-7** (2 x Xeon Gold 6130, Skylake-SP, 16 cores / 32 PU per socket, 2 NUMA nodes, one
22 MB L3 per socket, 6 channels of DDR4-2666 per socket, 187 GB; `oarsub -p grvingt --project cryptanalyse
-l walltime=3 "sleep 10800"`, job 6925944).  Commit `12fb293` (`omp_reboot` = `router-fat-slack`: the
`--prefetch` ring at default 8, the 16x pool slack, `block_points` 16384) on `bench-grvingt`, release build
(`-O3 -DNDEBUG -march=native -ggdb`), debian13, `module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0`.
Transparent huge pages `always`, `nr_hugepages` 0, so every shard came up on transparent 2 MB pages (the
report line confirmed full coverage, `thp_fault_fallback` stayed 0).  The kernel's NUMA balancing
(`kernel.numa_balancing` = 1 on this image) was scanning the shards (`numa_pte_updates` 159 M, 287 K hinting
faults after 20 runs) and was switched off after the first batch, section 8.

**The question:** what `--prefetch` buys on an Intel machine with a fraction of grdix's cores, where each
thread spends its time, what the limiting factor is, and whether a prefetch in FILL would help too.  The
run is session 15's, at this machine's size: `double_speck64_demo --n 40 --ram 160G --producers-per-node P
--dicts-per-node R`, one rank, `1 + P + R = 32` threads pinned one per core (the second hyperthread of every
core idle except in section 6), killed 30 s into the first PROBE (`build/session17/run17.sh`, the guarded
runner of session 15 with the roles' CPUs read from `/proc`).  160 GB is 21.5 G slots in R shards, fill 0.5:
a round inserts 1.0e10 = 2^33.2 preimages and probes 2^40; 110 rounds would make the search.  The
per-round collision rate is 1 % of probes (m = 40, 1.0e10 entries), against 5.5 % on grdix.

Every number is from one run unless said otherwise; the 8/23 configuration was run twice at
`--prefetch` 0 and 8 and repeated to 0.3 %.

## Summary

The dict thread and the producer each have a **fixed cost per point that does not depend on the split**
(2 % over P from 6 to 16), so the machine runs at min(P x producer rate, R x dict rate) and the team's
shape decides which role is the wall:

| per thread | PROBE, `--prefetch 0` | PROBE, `--prefetch 8` | FILL |
| --- | --- | --- | --- |
| **dict thread** | **15.1 M probes/s** (186 cycles, 82 instructions, IPC 0.44) | **22.5 M/s, +50 %** (124 cycles, 110 instructions, IPC 0.89) | 19.5-20 M inserts/s shipped; **26 M/s with a write-intent prefetch 16 ahead (variant CW), +33 %** |
| **producer, with the push** | 50 M points/s (48 cycles at the 2.4 GHz AVX-512 clock) | same | 47 M/s |
| **producer, never pushing** | 85 M/s (28 cycles) | same | 85 M/s |

So **`Router_Push` costs a producer 20 cycles a point here too** (grdix: 21), and the prefetch takes the
dict thread from one exposed DRAM miss per probe (73 % of its cycles on the slot load) to 124 cycles of
which about half are still spent waiting, section 4.  The best team with the prefetch is around 10 producers
and 21 dict threads (P x 50 = R x 22.5) and runs PROBE at about 0.47 G points/s; without the prefetch it is
around 7 / 24 at 0.35 G/s.  **At the split the prefetch does not change, its gain is what the receivers'
share of the wall allows: +14 % at 8/23 (the 8 producers cap the run at 400 M/s), +48 % at 12/19, +52 % at
16/15.**  DRAM is not the wall anywhere: PROBE moves 55-60 GB/s of reads out of a 256 GB/s peak.

## 1. The runs

| producers / dicts | FILL, 1.0e10 inserts | blocks held back in FILL | PROBE `--prefetch 0` | PROBE `--prefetch 8` | receivers discard (A) | producers never push (B) |
| --- | --- | --- | --- | --- | --- | --- |
| 6 / 25 | 35.4 s | 8.8 K | 296 M/s | 297 | -- | -- |
| 8 / 23 | 27.3 s | 61 K | **346** | **394** | 400 (FILL 26.7 s) | 680 (FILL 14.7 s) |
| 10 / 21 | 25.4 s | 450 K | -- | **461** | -- | -- |
| 11 / 20 | 26.8 s | 605 K | -- | 443 | -- | -- |
| 12 / 19 | 27.1 s | 608 K | 287 | **426** | -- | -- |
| 16 / 15 | 33.6 s | 608 K | 226 | 343 | 808 (FILL 14.9 s) | 1365 (FILL 7.4 s) |

PROBE is the live line's cumulative rate about 30 s in (the completed fraction of 2^40 over the elapsed
time agrees to 1 %).  The 8/23 sweep of `--prefetch`, in the interleaved order 0, 8, 4, 16, 2, 32, 0, 8,
64: **346, 395, 394, 395, 392, 394, 346, 393, 394 M/s** -- the knee is at 2 (as grdix's), 4 to 64 flat
within 1 %, and the whole plateau sits on the producers' ceiling of 400 M/s (A), so at this split it says
nothing about how far the dict thread can go.  At 12/19 and 16/15 the receivers are the wall (600 K blocks
held back for a full inbox, the producers 45 % of their time in the free-block wait) and the gain is the
dict thread's own: 15.1 -> 22.4 M/s and 15.1 -> 22.9 M/s per thread.

**Variants A and B** (`~/mitm-exp17` on the Nancy home, a scratch worktree at `12fb293` with `-DEXP_A`:
the dict thread counts the point and drops it; `-DEXP_B`: the producer computes, hashes and picks a
destination but folds the push into a sink) give the producers' ceilings: 50 M points/s per producer with
the push at P = 8 and 16 alike, 85 M/s without, in PROBE (running g) and FILL (running f) alike.  The
producers run at **2.38-2.39 GHz** under the AVX-512 licence where the dict threads run at 2.79 (perf's
cycles over 5 s, per role), so 20 ns a point with the push is 48 cycles, 11.8 ns without is 28, and the
push is 20 cycles -- session 15's 21 on grdix, at 6.8 ns there and 8.2 ns here.  `--benchmark` on this
machine was not run; session 15's caveat (dependent chains, 4x low) stands.

## 2. Where each thread spends its time (perf, 15 s of PROBE at 12/19, 499 Hz, every CPU)

The placer pins, so `-C` is an exact role filter (the omp order is service, R dict threads, P producers;
`run17.sh` writes each run's thread-to-CPU map).  Counters are `perf stat -C <role>` over 5 s of the same
PROBE, divided by the points the live line says the phase retired in those 5 s.

**Dict thread, `--prefetch 0`** (286 M/s = 15.1 M/s a thread, the wall): **186 cycles a point, 82
instructions, IPC 0.44**, 11.7 branches of which **1.12 mispredicted**, 3.2 L1D misses; a miss is
outstanding **86 % of the cycles** and 2.6 misses are in flight on average when one is
(`l1d_pend_miss.pending` / `pending_cycles`).  By source line: `dict.hpp:144`, the slot load of
`probe`'s run loop, **72.8 %**; `:184`, the counter after the probe, 8.9 % (skid); `:145` the tag compare
6.4 %; `murmur64` 1.9 %; `Router_Pop` 1.0 %; scalar Speck for the 1 % of probes that collide 0.8 %.  One
DRAM miss per point with nothing behind it, the same picture as grdix's 330 cycles, at this machine's
latency.

**Dict thread, `--prefetch 8`** (427 M/s = 22.5 M/s a thread, still the wall): **124 cycles a point, 110
instructions, IPC 0.89**, 13.2 branches, 1.10 mispredicted, 2.7 L1D misses; a miss outstanding **61 % of
the cycles**, 1.5 in flight when one is.  By line: `:144` the slot load **40.8 %**; **`:157`, the
`__builtin_prefetch` itself, 14.5 %**; `:184` 10.6 %; `Router_Pop` (`workers.hpp:241`) 5.5 %; the tag
compare 5.1 %; the ring's copy and swap (`:224-234`, `stl_pair.h`) 6.5 %; `murmur64` 4.3 %; scalar Speck
1.1 %.  The 28 instructions the ring adds are a quarter of the instruction count and cost about 30 cycles
at this IPC; the rest is that the slot load still waits.  Two candidate reasons, both measurable and
measured in section 7: the second-level TLB (1536 entries of 2 MB cover 3 GB of a 7 GB shard, so about
half the probes walk the page table, and a software prefetch that misses the TLB walks it too -- the 14.5 %
on the prefetch instruction is what that would look like), and prefetches arriving late.

**Producer at 12/19, `--prefetch 8`** (throttled: 12 x 50 = 600 M/s possible, 427 taken): 23 % in the
`pause` of the free-block wait (`tools.hpp:32`, `workers.hpp:86`), then the same profile as when it is
the wall.  **Producer as the wall** (8/23, `--prefetch 8`, 394 of its 400 M/s): **48.5 cycles a point, 98
instructions, IPC 2.03**, 1.3 L1D misses of which 0.36 are L2 hits, no L3 misses.  By line: vector Speck
(`double_speck64_problem.hpp:48-50, 58, 65` and the AVX-512 intrinsics) **48 %**; **`Router_Push`
(`workers.hpp:158-165`) 27 %** -- the point store `:161` 11.8 %, the count-word load `:160` 9.3 %, the rest
6 %; `murmur64` 7.9 %; the lane loop and destination pick (`producer.hpp:39-43`) 7 %; the f/g blend
(`:137`) 2.7 %.  The push's 27 % of the profile is 41 % of the producer's time by the intervention
(1 - 50/85), the same 3x-vs-profile gap as grdix's 11.5 % against 37 %, smaller here; the count word's
load being the second line of the push is what session 15's L1 hypothesis predicts (the count sits in
the last point of a 4 KB buffer, one line per destination).

**Service thread** (np 1, no peer): `Router_Progress` 24 %, Open MPI's `Testsome` and its mutexes
**about 40 %** (`ompi_request_default_test_some`, `pthread_mutex_lock/unlock`, `ompi_mtl_ofi_progress`),
`clock_gettime` 4.5 %.  A spinner at 22.6 K blocks/s against session 14's 2.0 M/s ceiling; off the
critical path, as on grdix.  (The MTL is `ofi` here: the Guix Open MPI picked libfabric over the
Omni-Path device even with one rank -- irrelevant at np 1, and the transport to force at np > 1 is
session 10's `--mca pml ucx`.)

## 3. DRAM traffic: not the wall, and twice the probe's line

`perf stat -a -e uncore_imc/cas_count_read/,uncore_imc/cas_count_write/` over 5 s of PROBE:

| run | reads | writes | per point |
| --- | --- | --- | --- |
| 8/23 `--prefetch 0`, 345 M/s | 60.4 GB/s | 9.2 GB/s | 175 B read, 27 B written |
| 8/23 `--prefetch 8`, 394 M/s | 57.7 GB/s | 10.5 GB/s | 146 B read, 27 B written |
| 12/19 `--prefetch 0`, 286 M/s | 49.8 GB/s | 7.7 GB/s | 174 B read, 27 B written |
| 12/19 `--prefetch 8`, 427 M/s | 62.8 GB/s | 11.5 GB/s | 147 B read, 27 B written |

A probe is one 64-byte line plus 16 bytes of block read; the machine reads **2.3-2.7 lines a point**.
The extra line is consistent with Skylake's L2 adjacent-line prefetcher pairing every random line with
its 128-byte buddy (hypothesis: not tested; MSR 0x1a4 bit 1 turns it off and would settle it in one run).
Either way 70 GB/s is a quarter of the 12 channels' 256 GB/s peak: the wall is latency and the
thread's own instructions, not bandwidth -- the opposite of grdix's FILL, where 192 dict threads did reach
the random-write ceiling (session 16).

## 4. FILL: a per-thread latency limit here, and a prefetch buys a third of it

FILL took **27.1-27.3 s at 8/23, 12/19 and 11/20 alike** (0.37 G inserts/s), 25.4 s at 10/21 (0.39), 33.6 s
at 16/15 (0.30), 35.4 s at 6/25 (0.28), for every `--prefetch` -- as it must, since `dict_round` inserts
straight from the pop in FILL and rings only the probes (the ring was removed from FILL after grdix's
sessions 15-16, where 192 dict threads sat on the machine's random-write wall).  But **this machine is
nowhere near a write wall**: at 6/25 and 8/23 the producers are the limit (47 M points/s each with the
push, few blocks held back), and at 12/19 and 16/15, where 600 K blocks are held back for a full inbox,
the receivers are, at **19.5 M inserts/s at 19 threads and 20 M/s at 15** -- the same rate per thread with
15 or 19 of them, a per-thread limit, with the dictionary's random writes at 0.3-0.4 G lines/s against a
memory system that reads 60 GB/s of random lines in PROBE without noticing.  (Inserts run faster than
unprefetched probes, 51 ns against 66, because FILL's runs are walked at an average load of 0.25 where
PROBE's are at 0.5.)

So the ring was put back into FILL in the scratch worktree, two ways: **C**, the shipped `dict.prefetch`
(read intent) before the insert too, and **CW**, `prefetchw` (write intent, `__builtin_prefetch(p, 1)`)
for inserts.  At 12/19, `STOP=fill`:

| FILL at 12/19 | 1.0e10 inserts | inserts/s per dict thread | blocks held back |
| --- | --- | --- | --- |
| shipped, any `--prefetch` | 27.1 s | 19.5 M | 608 K |
| C, `--prefetch 0` (the same loop) | 28.5 s | 18.5 M | 598 K |
| C, read intent, 8 ahead | **21.7 s** | 24.3 M | 526 K |
| C, read intent, 16 ahead | 21.3 s | 24.7 M | 373 K |
| CW, write intent, 8 ahead | 20.5 s | 25.7 M | 258 K |
| **CW, write intent, 16 ahead** | **20.3 s** | **26.0 M** | 219 K |

**A prefetch in FILL is worth 25-33 % on this machine**, write intent beating read intent by 5 % (the
line arrives in the right state for the store, no second transition), 16 ahead a hair over 8; and the
blocks held back fall 3x, so at 20 s the receivers are close to the 12 producers' own FILL ceiling of
about 18 s (A: 47 M/s x 12).  C's own `--prefetch 0` at 28.5 s against the shipped binary's 27.1 is the
run-to-run scatter of this session's worst pair, 5 %, or the variant's extra `D == 0` test; the gain
is measured against both.  The micro-benchmark of section 5 says the same thing without the Router.

This is the opposite verdict from grdix (sessions 15-16: FILL unmoved by any prefetch), and both are
right: there the wall was the machine's, 192 threads x 10.8 M inserts/s = 2.07 G random dirty lines/s
against a 3.35 G/s ceiling; here it is the thread's, 19 x 19.5 M = 0.37 G/s against a DDR4 system that
has not been asked for a tenth of what grdix's was.  **The FILL prefetch should be shipped and be the
same knob** (`--prefetch`, write intent for inserts): it costs nothing where the wall is the memory and
a third of FILL where it is the thread.

## 5. The dict thread alone: `fillbench`

`build/session17/fillbench.c`: T OpenMP threads pinned one per core (`OMP_PLACES=cores
OMP_PROC_BIND=spread`), each owning a private 4 GB shard on transparent 2 MB pages, fed a synthetic
stream of uniformly random keys through `murmur64`; the FILL inserts 0.5 x n_slots entries with
`DirectDict::insert`'s loop, the PROBE probes as many random keys against the filled shard with
`probe`'s, each through the same ring of D points as `dict_round`, each thread timing its own phase.  No
Router, no producers, no block stream: what the dict thread's loop costs by itself.

Points/s per thread, the phase time being the slowest thread's; T = 23 is the dict count of the 8/23
team, 31 every worker core, 8 and 1 the uncontended thread (`batch2.out`, 09:14-09:27):

| D ahead | FILL, T = 23 | FILL, T = 31 | FILL, T = 1 | PROBE, T = 23 | PROBE, T = 31 | PROBE, T = 1 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 25.0 M/s (39.9 ns) | 24.1 | 29.5 | **15.2 M/s (65.7 ns)** | 14.9 | 17.0 |
| 2 | 25.1 | 22.8 | 33.6 | 20.8 | 20.1 | 24.0 |
| 4 | 29.0 | 25.4 | 40.5 | 24.2 | 23.5 | 28.3 |
| 8 | 34.1 (+36 %) | 29.2 | 49.7 | **27.1 (+78 %)** | 26.1 | 31.9 |
| 16 | 38.2 (+53 %) | 33.2 | 58.0 | 27.9 | 26.8 | 33.5 |
| 32 | 40.1 (+60 %) | 35.2 | 61.5 | 27.5 | 26.4 | 33.2 |
| 8, write intent | 35.7 | -- | -- | (27.1) | -- | -- |
| 16, write intent | 38.9 | -- | -- | (27.9) | -- | -- |

The unprefetched probe alone is **15.2 M/s a thread, the engine's 15.1 to the decimal**: the dict thread
in PROBE is its probe loop and nothing else.  Prefetched, the loop alone reaches 27.1 M/s where the
engine's dict thread reaches 22.5: the 7 ns between them are `Router_Pop`, the block stream and the
collision check, 17 % of the prefetched thread's 44 ns.  The knee is the engine's (2 buys half, 8 nearly
all; PROBE is flat from 8 on, FILL keeps climbing to 32 because an insert is the shorter operation).
**FILL alone gains 36-60 %**, write intent 2-5 % more, confirming section 4 from outside the tree; and at
T = 31, D = 32 the machine writes 1.09 G lines/s, which is **the random-write ceiling `randread2` measures
(1.05-1.07 G lines/s at 32 threads, any number of streams)** -- 31 prefetched dict threads would sit on
this machine's write wall in FILL just as 192 unprefetched ones did on grdix, at a third of grdix's
line rate.  With 19-21 dict threads the engine stays under it (0.49 G/s in the CW run).

The load matters more than the prefetch: at fill 0.25 a prefetched probe costs **20 ns (50 M/s a
thread)**, at 0.5 it costs 37, at 0.75 it costs 76 -- the run walked behind the prefetched first line,
crossing into a second line the prefetch did not fetch, is what the prefetched thread pays for.  A layout
that keeps a key's run inside one 64-byte line (8-slot buckets, or the 8-slot compare of session 16) would
give the fill-0.5 dictionary the fill-0.25 probe.

**The machine's random-access ceilings (`randread2`, 32 threads x 3 GB, 2 MB pages):** reads with
independent addresses **1.65-1.76 G lines/s** (113 GB/s of lines, 44 % of the 12 channels' peak) from one
stream a thread up -- a Skylake core overlaps them by itself at 52-55 M lines/s, 18-19 ns each, which is
also the per-thread limit at 23 threads (1.28 G/s); dependent addresses at one miss in flight **0.33 G/s,
98 ns**, at 8 in flight 1.73 G/s; **writes 1.05-1.07 G lines/s** whatever the number of streams.  The
unprefetched dict thread (15 M/s) does better than the one-miss chase (10 M/s) because the core overlaps
the tail of a probe with the head of the next, 2.6 misses in flight (section 2).  And the best PROBE of
this session, 613 M/s with two threads per core, at 2.3-2.7 lines a point is **1.4-1.6 G lines/s, 80-90 %
of the random-read ceiling**: once the prefetch and the second hyperthread are both in, DRAM becomes the
next wall on this machine, and the adjacent-line prefetcher's doubling of the traffic (section 3,
hypothesis) becomes worth turning off.

## 6. The team's shape, and the second hyperthread

With the per-thread rates of section 1 the best one-thread-per-core team is where P x 50 = R x 22.5 with
P + R = 31: **10 producers and 21 dict threads**, measured **461 M points/s** in PROBE (11/20: 443, 12/19:
426, 8/23: 394 at the producers' ceiling).  Without the prefetch the balance is near 7/24 and 350 M/s; the
prefetch is worth about **+32 % at the best split of each**, and moves the split by three cores.

The 64 PUs were then all used, two threads per core, the placer pairing them by its emptiest-core rule
(one thread per core first, then the second PU): at **24 producers + 39 dict threads** every producer
shares its core with a dict thread (24 dict+prod cores, 7 dict+dict, 1 service+dict), no two producers
share the AVX-512 units of one core, and the dict threads' misses overlap with the sibling's Speck:

| team (P/R) | threads per core | FILL | PROBE `--prefetch 0` | PROBE `--prefetch 8` |
| --- | --- | --- | --- | --- |
| 10 / 21 | 1 | 25.4 s | -- | 461 M/s |
| 16 / 47 | 2 | **20.9 s** | 468 | 562 |
| 24 / 39 | 2 | 22.6 s | 405 | **613 M/s** |

**Two threads per core is worth +33 % in PROBE and +18 % in FILL over the best one-per-core team**, with
the shipped binary: a hyperthreaded dict thread runs at about 12-16 M probes/s instead of 22.5, but there
are 39 of them, and the producer beside each keeps 25 M points/s of the 50.  That is the same lever as
the prefetch, misses in flight per core, bought with the other PU instead of a ring; and the two stack
(613 against 405).  The grdix reference team (63/192, one per core, 512 PUs idle) has not tried this.
Not measured: 20/43, 28/35 and the SMT team with the FILL prefetch.  Section 5's ceilings say the SMT team is already at 80-90 % of the machine's random-read rate.

## 7. What the prefetched dict thread still waits for: the TLB and the run's tail

`perf stat` on the 19 dict CPUs, 12/19, 5 s of PROBE at 428 M/s (2.16 G points), `--prefetch 8`:
`dtlb_load_misses.miss_causes_a_walk` **0.63 a point**, all of them 2 MB walks, `walk_active` **16 % of the
cycles = 20 cycles a point**, `stlb_hit` 0.76 a point; `sw_prefetch_access.t0` 1.16 a point and
`load_hit_pre.sw_pf` **0.04 a point**: the prefetches do arrive in time, only 4 % of loads find their line
still in flight.  `cycle_activity`: **48 % of cycles stalled**, 40 % on memory, **23 % on an L3 miss** (28
cycles a point), 37 % on an L1D miss; `offcore_requests_outstanding` 5.2 data reads in flight when any,
95 % of cycles with one.  Without the prefetch: 0.66 walks a point, `walk_active` 10 % (18 cycles), and
the rest is the exposed miss.

So the 124 cycles are about 60 of instructions (110 at the IPC of the unstalled cycles), **20 of page
walks** -- the 7 GB shard is 3584 pages of 2 MB against a 1536-entry STLB, so two probes in three walk,
into a 28 KB PMD table that sits in L2 -- and **28 of L3-miss stall that the prefetch does not cover**:
the second line of a run that crosses one (a run of 2.5 slots at load 0.5 crosses a 64-byte line 1 time
in 4), the block stream's line every 4 points, and the prefetch that was dropped or came too late (4 %).
The levers, in order: fewer instructions (carry the hash from the prefetch to the probe, retire in place,
compare a line of 8 slots at once -- session 16's list), 1 GB pages for the shards (the walks: needs
`hugepagesz=1G` at boot, untested), prefetching the *next* line of the run too when the first slot is
near a line's end.  1 GB pages and the 8-slot compare are the two that Skylake would feel most.

## 8. Method

- **NUMA balancing was on** (`kernel.numa_balancing` = 1 on the debian13 image) and scanning the shards:
  one run in 22 filled in 62 s instead of 27, the slowdown starting 10 s in and growing, with THP coverage
  intact; its PROBE matched its twins.  Off for the rest of the session (`sudo-g5k sysctl -w
  kernel.numa_balancing=0`); the Router pins and first-touches, so balancing can only move pages away
  from their thread.  Turn it off before measuring on any node.
- 160 GB of 187: `MemAvailable` came back above 178 GB within a second of every `kill -9` here, unlike
  grdix's 960 GB (session 15); the guard stays.
- A `nohup ... &` inside an `ssh` command holds the ssh session until the job ends, even with every
  stream redirected and `setsid`: launch the batch, then read its `.out` file from the front-end (the home
  is NFS), never through the launching connection.
- The frontend has the same `perf` as the nodes and sees the same files: `perf report` of a node's
  `perf.data` runs there without touching the node under measurement.
- Runs at a producer-bound split cannot measure the receiver (the 8/23 sweep's flat plateau is the
  producers' ceiling); pick the split so that the role under study is the wall, and check it is (blocks
  held back, the other role's `pause` share).

## What to change

1. **Ship the FILL prefetch, write intent, on the same `--prefetch`**: +33 % of FILL here (27.1 -> 20.3 s),
   0 % on grdix where the wall is the memory's, so nothing to lose.  Section 4.
2. **The team's shape is the first-order knob on a 32-core machine**: 10/21 (461 M/s) against 8/23 (394)
   and 16/15 (343) with the prefetch; the per-thread rates (50 with the push, 22.5 prefetched, 15 not)
   predict it.  Consider deriving the default split from them instead of asking for both counts.
3. **Use the second hyperthread**: 24/39 two per core gives 613 M/s, +33 % over the best one-per-core
   team, the placer already pairs a producer with a dict thread on each core.  Try the SMT team on grdix
   (63/192 leaves 256 PUs idle) and the SMT team with the FILL prefetch here.
4. **The push is 20 cycles on Skylake too** (41 % of a producer's time when it is the wall), the count word
   the second-hottest line of the push: session 15's compact-counts test is still the one to run.
5. **The prefetched dict thread's 124 cycles**: 60 instructions, 20 page walk, 28 L3 stall.  Instructions
   first (session 16's list), then 1 GB pages, then the run's second line.  Section 7.
6. `kernel.numa_balancing = 0` in the run scripts; and the `--prefetch` default of 8 is right here too
   (knee at 2, flat 4-64).

## Reproducing

```bash
oarsub -p grvingt --project cryptanalyse -l walltime=3 "sleep 10800"    # then ssh grvingt-N.nancy.g5k
module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0
cmake -S . -B build -DCMAKE_BUILD_TYPE=release \
    -DHWLOC_INCLUDE_DIR=/gnu/store/hl3isj962ghsvvxig6ls84mw54nfc917-hwloc-2.13.0-lib/include \
    -DHWLOC_LIBRARY=/gnu/store/hl3isj962ghsvvxig6ls84mw54nfc917-hwloc-2.13.0-lib/lib/libhwloc.so
make -C build -j 32 double_speck64_demo
sudo-g5k sysctl -w kernel.perf_event_paranoid=-1 kernel.kptr_restrict=0 kernel.numa_balancing=0
# one run: P producers, R dicts, killed 30 s into PROBE (STOP=fill: 3 s after the FILL report); PERF=1 STAT=1 profile it
EXTRA="--prefetch 8" TAG=p10r21-pf8 P=10 R=21 bash -l build/session17/run17.sh
build/session17/summ17.sh p10r21-pf8 ...             # FILL s, blocks held back, PROBE rate, rate from the completed fraction
build/session17/report17.sh p12r19-pf8p 19           # per-role profiles by symbol and source line (runs on the front-end too)
# the variants: ~/mitm-exp17 (A: receivers discard, B: producers never push, C / CW: the ring in FILL, read / write intent)
cmake -S ~/mitm-exp17 -B ~/mitm-exp17/build-CW -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS=-DEXP_CW ...
# the batches of this session: batch1.sh (ceilings, sweep, splits, profiles), batch3.sh (FILL variants, 12/19 profiles,
# 10/21, 11/20, SMT), batch4.sh (TLB and stall counters), batch2.sh (fillbench + randread2, the dict thread alone)
```

`build/session17/` in the grvingt worktree keeps every log, `.tids` map, `perf.data`, stat log and script.

# Session 19 -- 8 grvingt hosts on Omni-Path: the direct engine is bound by the service thread's SDMA page work, and the hfi1 pin cache buys 19 %

2026-09-15, **grvingt-1, 2, 4, 6, 7, 8, 13, 24** (2 x Xeon Gold 6130, 32 cores / 64 PU, 187 GB, Omni-Path 100 Gb/s),
`oarsub -t deploy -p grvingt --project cryptanalyse -l nodes=8,walltime=3:00` (job 6926976), **kadeploy of
`debian13-big` with `iommu.passthrough=1` in its kernel parameters** (`~/debian13-big-iommupt.yaml`, the durable
form of sessions 12-13's kexec), kernel 6.12.107, `module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0 gcc-toolchain`
(Guix Open MPI 4.1.6, PSM2), `kernel.numa_balancing = 0`.  Commit `4d4d847` on `bench-grvingt` = `fe7ae49` of
`omp_reboot` plus session 17's FINDINGS.  Release build on the node.  Scripts in `build/session19/`.

**The question** (the user): the double_speck64 benchmark over 8 hosts with 32 threads per node, and whether the
same rate is reached with fewer threads because the network is the bottleneck.

**The run**: `double_speck64_demo --n 38 --ram 40G --prefetch 8 --nrounds 1 --producers-per-node P
--dicts-per-node R`, one rank per host, the Router pinning its team on the first hardware thread of each core,
`mpirun --hostfile nodes --map-by ppr:1:node --bind-to none -x PATH -x HFI_NO_CPUAFFINITY=1 --mca btl ^openib
--mca pml cm --mca mtl psm2` from grvingt-1 in a login shell (`run19.sh`).  `--nrounds 1` ends the run after round
0's PROBE with the exact round report (the live line on several nodes is rank 0's share scaled by 8, and the
FILL-to-PROBE percentages it prints are that approximation).  FILL is 2^34.2 inserts a round, PROBE 2^38 probes;
the engine's `node-->` figure is the bytes each node sends per second.

## Summary

| team per node (P / R) | threads | PROBE, G points/s, 8 hosts | per node | wire, GB/s per node |
| --- | --- | --- | --- | --- |
| 10 / 21, hfi1 pin cache 256 MB (the default) | 32 | 2.22 | 0.28 | 3.9 |
| **10 / 21, pin cache 4 GB** (two runs) | 32 | **2.64, 2.40** | 0.33, 0.30 | 4.6, 4.2 |
| 10 / 21, pin cache 4 GB, 1 MB blocks | 32 | 2.36 | 0.29 | 4.1 |
| 10 / 21, pin cache 4 GB, `--n-recv 128 --credit 8` | 32 | 2.53 | 0.32 | 4.4 |
| 8 / 15, pin cache 4 GB | 24 | 1.97 | 0.25 | 3.5 |
| 7 / 12, pin cache 4 GB | 20 | 1.77 | 0.22 | 3.1 |
| 6 / 9, pin cache 4 GB | 16 | 1.31 | 0.16 | 2.3 |
| 5 / 10, pin cache 4 GB | 16 | 1.41 | 0.18 | 2.5 |
| 4 / 7, pin cache 4 GB | 12 | 1.09 | 0.14 | 1.9 |

**The 8-host engine runs at 0.33 G points/s per node, 72 % of the same team's 0.46 on one node (session 17), with
the wire at 4.6 of its 11.6 GB/s.**  It is neither the producers nor the dict threads nor the fabric: both worker
roles spend a fifth of their samples in `cpu_relax`, and the one thread that moves the blocks over MPI -- the
service thread -- spends most of its time in the kernel's hfi1 SDMA path pinning the pages of every message.
The driver's pin cache (`/sys/module/hfi1/parameters/cache_size`, 256 MB) is smaller than the Router's pool
(1.9 GB), so the cache thrashed; **4 GB is +19 % on PROBE and +32 % on FILL**, and the receive side never runs
out of free blocks again.  Bigger blocks and more in-flight receives and credit do nothing or hurt: the service's
cost is per byte, in the kernel.  **Fewer threads do not reach the maximum**: with the delivery fixed, PROBE scales with the dict threads (16-18 M
probes/s each, whatever the team), 24 threads lose a quarter and 16 threads half; only FILL saturates early, at
2.4-2.5 G/s from 20 threads up, because FILL is where the service thread's per-byte cost binds.  **Identical runs
scatter 9 %** (10 / 21 with the cache: 2.64 then 2.40), so the cache's PROBE gain is +8 to +19 % and the knobs'
differences below are inside the noise; the FILL gain (+32 %, 2.50 both times) and the vanished "receives left
unposted" are the robust part.

## 1. Getting 8 nodes to run at all

Four things, each of which cost a launch:

1. **A deploy job's nodes are unreachable until deployed** (ssh as the user: "not all its CPU cores are assigned
   to the job, use oarsh"; `oarsh`: no cpuset; root: no key), so there is no kexec "out of the box"; and once
   deployed there is no need for one: `kadeploy3 -a ~/debian13-big-iommupt.yaml -f nodes -k` boots the kernel
   with `iommu.passthrough=1` directly (15 min for 8 nodes, two of which needed kadeploy's own hard reboot).
   `post19.sh` then checks `Default domain type: Passthrough` and does the OPA rename on every node; `flood19.sh`
   measured 11.6-11.7 GB/s over PSM2 between two of them (a cold first run 3.8), Open MPI's default choice
   included.
2. **libpsm2 pins the whole process to one CPU at `MPI_Init`.**  Every rank came up with a 1-CPU affinity mask
   whatever `--bind-to none` said (a `bash -c taskset` under the same mpirun shows 0-63), and the Router refused
   with "rank 0 has 1 CPUs in its affinity mask for 32 threads".  It happens under `--mca pml ucx` too (Open MPI
   probes the psm2 MTL anyway).  **`-x HFI_NO_CPUAFFINITY=1`** on the mpirun line fixes it.  It never showed
   before because PSM2 never worked before (sessions 10, 17).
3. The deployed image has no `br0`: `--mca oob_tcp_if_include br0` aborts the launch; drop it.  Redeployed nodes
   change host keys: `ssh-keygen -R` them, and launch with `--mca plm_rsh_agent "ssh -o StrictHostKeyChecking=no"`.
4. **A problem size that makes the dictionary a large share of the range inflates the collision check.**  The
   first sweep ran `--n 38 --ram 160G`: 2^36.2 entries in a 38-bit range means **29 % of probes hit a genuine
   collision** and pay `is_good_pair`'s scalar Speck (key schedule + encryption), against 5.5 % in the grdix
   reference (`--n 40`, 960G).  The dict threads then run 290 cycles a probe (35-50 % of their samples in Speck),
   the whole machine is receiver-bound and every number is about Speck.  Discarded (`ram160/`); `--ram 40G`
   brings the share to 7 % and FILL to 10 s.  The rule: `fill * w / 2^m` is the collision rate a probe pays for,
   keep it where the reference has it.

## 2. The team sweep (pin cache at its 256 MB default)

| P / R | threads | FILL s | FILL G/s | PROBE s | PROBE G/s | wire in PROBE, GB/s per node | blocks held back |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 / 4 | 8 | 23.1 | 0.87 | 527.8 | 0.52 | 0.9 | 1.0 M |
| 4 / 7 | 12 | 12.8 | 1.56 | 248.2 | 1.11 | 1.9 | 0.4 M |
| 5 / 10 | 16 | 11.0 | 1.82 | 195.1 | 1.41 | 2.5 | 1.1 M |
| 8 / 15 | 24 | 10.6 | 1.89 | 139.7 | 1.97 | 3.4 | 1.1 M |
| 8 / 23 | 32 | 10.6 | 1.89 | 125.7 | 2.19 | 3.8 | 1.1 M |
| **10 / 21** | 32 | 10.6 | 1.89 | 123.8 | **2.22** | 3.9 | 1.1 M |
| 12 / 19 | 32 | 10.6 | 1.89 | 130.4 | 2.11 | 3.7 | 1.1 M |
| 14 / 17 | 32 | 10.4 | 1.92 | 138.4 | 1.99 | 3.5 | 1.2 M |

The 32-thread teams are flat at 2.0-2.2 G/s whatever the split, the 24-thread one is within 12 % of them, and
below 16 threads the rate falls with the producers.  FILL is 1.9 G/s for every team of 24 threads or more.  A
plateau that the split does not move and that a smaller team nearly reaches is not a worker-side limit.

## 3. Where the time goes (perf on grvingt-2, 10 / 21, steady PROBE, pin cache 256 MB)

Per-role counters over 5 s (`perf19.sh`, root on the node, the roles from the pinned threads' CPUs):

| role | clock | instructions / point | cycles / point | in `cpu_relax` (`tools.hpp:32`) |
| --- | --- | --- | --- | --- |
| dict thread, 21 per node, 13.2 M probes/s each | 2.80 GHz | 177 | 212 | **19 %** of samples |
| producer, 10 per node, 27.8 M points/s each | 2.40 GHz (AVX-512) | 115 | 86 | **24 %** |
| service | 2.79 GHz | | | -- |

The dict thread probes at 177 instructions a point (150 in session 18, plus the 7 % of Speck), so it is not
starved outright, but a fifth of its samples are the pause of an empty inbox; the producer's 115 instructions are
its 100 of work plus the pause of `refill` waiting for a free block (`workers.hpp:86`, and the free ring's CAS at
`:65`).  Both roles wait on the same thing.  Every dict thread runs at the same rate (2.45-2.57 G instructions/s
per CPU): no imbalance between receivers.

**The service thread** (one per node, all the MPI) by symbol, 10 s of PROBE:

| share | where |
| --- | --- |
| 14 % | `Router_Progress` (its own loop) |
| 9.7 % | `__mmu_int_rb_subtree_search` (kernel: the hfi1 pinned-pages tree, looked up per SDMA packet) |
| 7.8 % | `hfi1_add_pages_to_sdma_packet` (kernel: the pages of a message into SDMA descriptors) |
| 5.5 % + 2.2 % | `ompi_request_default_test_some`, `PMPI_Testsome` |
| 3.3 % | `user_sdma_send_pkts` (kernel) |
| 3.7 % + 2.1 % + ... | `__list_del_entry_valid_or_report`, `memset_orig`, the syscall entry/return pairs |
| 1.5 % | `psm2_mq_ipeek2` |

and in the first, 160G profile also `unpin_user_pages`, `gup_fast_fallback`, `hfi1_mmu_rb_insert`: the pin
cache **evicting and re-pinning**.  PSM2 sends a 256 KB message by SDMA from the user's pages; hfi1 pins them and
keeps them in an MMU-notifier rb-tree whose size is the module parameter `cache_size` -- **256 MB**, while the
Router's pool is 7111 blocks = **1.9 GB** at 10 / 21.  Every message therefore re-pins its 64 pages and evicts
someone else's.  The Router counters agree: 1.1 M blocks held back for a full destination, and in PROBE
**1.5 M "receives left unposted for want of a free block"** -- the receiving service had no block to post a
receive into, the sending side's credits stalled, its sealed blocks parked, and its producers waited in
`refill`.  The user's reading ("the pool is exhausted because the network does not keep up with the senders, so
the blocks park") is the mechanism; what does not keep up is the one thread feeding the network, in the kernel.

## 4. The intervention: `cache_size` 256 MB -> 4 GB (`setcache19.sh`, root, writable at runtime)

Same team, same problem, one parameter:

| | FILL s | FILL G/s | PROBE s | PROBE G/s | wire, GB/s per node | receives left unposted in PROBE |
| --- | --- | --- | --- | --- | --- | --- |
| 256 MB | 10.6 | 1.89 | 123.8 | 2.22 | 3.9 | 1.5 M |
| **4 GB** | **8.0** | **2.50** | **104.2** | **2.64** | **4.6** | **0** |

+19 % on PROBE (+8 % on the repeat run: 2.40, the fabric's 9 % run-to-run scatter), +32 % on FILL (2.50 both times), and the receive side never runs dry again.  The dict threads' pause share fell
from 19 % to 9 % (`perf19.sh` during the 4 GB run), and the service thread's profile still leads with
`__mmu_int_rb_subtree_search` (13 %) and `hfi1_add_pages_to_sdma_packet` (10 %): the eviction is gone, the
per-packet lookup is not -- that is the cost of the SDMA path on a hit, and it scales with bytes.

**The Router's own knobs, on top of the 4 GB cache**: `--block 65536` (1 MB messages) 2.36 and FILL 2.1 (FILL
**loses** 16 %, a robust difference), `--n-recv 128 --credit 8` 2.53, both together 2.32 -- all inside the 9 %
run-to-run scatter on PROBE except the 1 MB blocks' FILL.  Fewer, bigger messages do not help because the
kernel's work is per page, and more in flight does not help because the service is not waiting on the fabric.
The 256 KB default stands, as session 10 found on this fabric.

## 5. Fewer threads, with the fix

Same problem, pin cache at 4 GB, the team shrunk at the same 1 : 2 ratio (`batch19.sh`, `TAGSUF=c`):

| P / R | threads | FILL s | FILL G/s | PROBE s | PROBE G/s | per dict thread, M probes/s | wire in PROBE, GB/s per node |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 / 21 | 32 | 8.0, 8.1 | 2.50, 2.50 | 104.2, 114.7 | **2.64, 2.40** | 15.7, 14.3 | 4.6, 4.2 |
| 8 / 15 | 24 | 8.5 | 2.40 | 139.4 | 1.97 | 16.4 | 3.5 |
| 7 / 12 | 20 | 8.2 | 2.40 | 155.7 | 1.77 | 18.4 | 3.1 |
| 6 / 9 | 16 | 10.3 | 2.00 | 209.5 | 1.31 | 18.2 | 2.3 |
| 5 / 10 | 16 | 10.6 | 1.90 | 195.3 | 1.41 | 17.6 | 2.5 |
| 4 / 7 | 12 | 12.5 | 1.60 | 252.4 | 1.09 | 19.5 | 1.9 |

**PROBE**: the dict thread runs at 14-19 M probes/s whatever the team (124 cycles a probe in session 17 with 5.5 %
collisions and local blocks, about 175 here with 7 % and blocks off the wire), so the aggregate is the dict-thread
count times that, and 24 threads are 25 % short of 32.  Below 24 threads nothing changed with the pin cache
(8 / 15: 1.97 before and after; 5 / 10: 1.41 both) -- those teams never touched the service's limit.  **FILL**: 2.4-2.5
G/s from 20 threads up, the same for 7 / 12 as for 10 / 21: that is the service thread's ceiling in the phase that
sends every point off the node, and the extra dict threads do nothing for it.  So, to the user's question: on this
fabric the network is not what bounds the engine -- the wire idles at 4.6 GB/s -- and the 32-thread team is
needed for PROBE; the one thread whose time is scarce is the service, and it binds FILL first.

## What to change

1. **On grvingt, boot `iommu.passthrough=1` (deploy `debian13-big-iommupt`), set `hfi1 cache_size` to cover the
   Router's pool (4096 MB), rename the OPA device, and launch with `-x HFI_NO_CPUAFFINITY=1`** -- four per-boot
   or per-launch settings, none of them code, worth respectively "works at all", +19 %, "PSM2 at all" and "the
   Router pins at all".  `post19.sh` and `setcache19.sh` do the node side; the durable form of the cache size is
   a `modprobe.d` option in the environment's postinstall.
2. **The multi-node engine is bound by the service thread's kernel time per byte sent**, 4.6 GB/s per node with
   the cache fixed against 11.6 on the wire.  Levers, none measured yet: a second service thread per node (the
   Router has one, by design), PSM2's eager/PIO path for the block size (`PSM2_MQ_RNDV_HFI_THRESH`,
   `PSM2_MQ_EAGER_SDMA_SZ`: a CPU copy instead of pinning), huge pages under the pool (fewer pages per message
   to look up, if hfi1 pins compound pages as one), or -- the session 8 remark -- keeping points on the node.
3. **Problem-size hygiene for multi-node runs**: keep `fill * w / 2^m` near the reference's 5 %, or the measurement
   is of `is_good_pair`.  `--n 38 --ram 40G` on 8 grvingt nodes is the setting that does.
4. **Keep the 32-thread team on grvingt** (10 / 21, one thread per core): PROBE scales with the dict threads to the
   end; FILL is service-bound from 20 threads up, so it is the phase where a cheaper send path pays first.  The
   8-host reference: `--n 38 --ram 40G --producers-per-node 10 --dicts-per-node 21 --prefetch 8`, PROBE 2.4-2.6 G
   points/s, FILL 2.5 G inserts/s, pin cache 4 GB; three interleaved repeats before believing a knob (9 % scatter).

## Reproducing

```bash
# frontend: reserve, deploy, verify, cache
oarsub -t deploy -p grvingt --project cryptanalyse -l nodes=8,walltime=3:00 "sleep 10800"
kadeploy3 -a ~/debian13-big-iommupt.yaml -f nodes.txt -k
bash build/session19/post19.sh          # cmdline, Passthrough, OPA rename, link state, per node
bash build/session19/setcache19.sh 4096 # hfi1 pin cache, per node
bash build/session19/flood19.sh         # 11.6 GB/s over PSM2 between the first two nodes
# grvingt-1, login shell: the sweep, the knobs, the profile
RAM=40G CONFIGS="10/21 5/10 12/19 8/15 8/23 4/7 14/17 3/4" bash -l build/session19/batch19.sh
P=10 R=21 N=38 RAM=40G KNOBS="cache4g|;cache4g_blk64k|--block 65536" bash -l build/session19/batch19k.sh
python3 build/session19/summ19d.py d10_21 k10_21_cache4g ...
ssh root@grvingt-2 bash build/session19/perf19.sh   # per-role counters and profile of the newest run, on the node
```

# Session 24 -- 16 grvingt hosts on Omni-Path: the wire at 12.2 GB/s, the engine at 0.30 G points/s a host with one rank, 0.45 with a rank per socket

2026-09-16, **grvingt-1, 2, 3, 4, 5, 7, 8, 25, 27, 28, 30, 31, 32, 34, 35, 52** (2 x Xeon Gold 6130, 32 cores / 64 PU, 187 GB,
Omni-Path 100 Gb/s), `oarsub -t deploy -p grvingt --project cryptanalyse -l nodes=16,walltime=3:00` (job 6928266),
**kadeploy of `debian13-big-iommupt`** (`~/debian13-big-iommupt.yaml`: `debian13-big` with `iommu.passthrough=1` in its
kernel parameters, sessions 12-13's finding in its durable form), kernel 6.12.107, `module load openmpi/4.1.6
hwloc/2.13.0 ucx/1.20.0 gcc-toolchain` (Guix Open MPI 4.1.6, PSM2), hfi1 `cache_size` 4096 MB and
`kernel.numa_balancing = 0` on every node.  Commit `922ca44` on `bench-grvingt` = `0bff8a2` of `omp_reboot` (the placer
keeps the service thread alone on its core) plus the FINDINGS of sessions 17 and 19.  Release build on grvingt-1.
Scripts in `build/session24/`.

**The questions** (the user): the 16-host double_speck64 benchmark on the IOMMU-free kernel, after a check that the
wire is near its 12.5 GB/s; and two ranks a node, one pinned to each socket.

**The run**: `double_speck64_demo --n 39 --ram 40G --prefetch 8 --nrounds 1`, 16 hosts: 80 G slots (2^36.22), 40 G
entries a round, `fill * w / 2^m` = 7.3 % (the 8-host setting of session 19 at twice the range, same rate), 14 rounds
for an exhaustive search, one run.  FILL is 2^35.22 inserts, PROBE 2^39 probes; the rate is the exact round report's
evaluations over its time, `summ24.py`.  One rank a host is `--map-by ppr:1:node --bind-to none`; two ranks a host is
**`--map-by ppr:1:socket --bind-to socket`**, 20G a rank, the Router pinning inside the socket's 16 cores (32 hardware
threads) with the service alone on its core: 15 cores for the workers, so 5/10 one a core or 10/20 two a core.  Four
ranks a host is `--map-by ppr:2:socket:pe=8 --bind-to core`, 10G a rank, 7 cores for the workers.  `run24.sh` takes
P, R, N, RAM, PPN, EXTRA; `batch24.sh` a list of configs; the batches were chained (`chain24*.sh`).

## Summary

| ranks a host | team a rank (P / R) | knobs | FILL, G inserts/s | PROBE, G points/s | PROBE a host | wire a host in PROBE, GB/s |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 10 / 21 (the session 19 reference) | -- | 4.00, 4.08, 4.13 | 4.48, 4.52, 4.84 | 0.28-0.30 | 4.2-4.5 |
| 1 | 10 / 21 | `--n-recv 64` / `128` | 4.71 / 4.55 | 4.82 / 4.50 | 0.30 / 0.28 | 4.5 / 4.2 |
| 1 | 24 / 38 (SMT) | -- | 3.85 | 4.37 | 0.27 | 4.1 |
| 2 | 5 / 10 | -- | 6.45, 6.45 | 3.51, 3.53 | 0.22 | 3.4 |
| 2 | 5 / 10 | `--n-recv 128` / `128` / `256` | 6.45 / 6.16 / 6.56 | 3.52 / 3.58 / 3.52 | 0.22 | 3.4 |
| 2 | 6 / 9 | -- / `--n-recv 128` | 5.07 / 5.07 | 3.75 / 3.87 | 0.23-0.24 | 3.6-3.8 |
| 2 | 10 / 20 (SMT) | -- | 4.71, 4.71 | 6.42, 6.29 | 0.39-0.40 | 6.0-6.2 |
| 2 | 10 / 20 (SMT) | `--n-recv 128` | 6.35, 6.56, 6.67 | 6.14, 6.85, 6.33 | 0.38-0.43 | 6.0-6.6 |
| 2 | 10 / 20 (SMT) | `--n-recv 128 --credit 8` | 4.88 | 6.58 | 0.41 | 6.4 |
| **2** | **8 / 22 (SMT)** | **`--n-recv 128`** | **6.45, 6.25, 6.35** | **7.07, 7.23, 7.14** | **0.44-0.45** | **7.0** |
| 2 | 7 / 23 (SMT) | `--n-recv 128` | 5.80 | 6.60 | 0.41 | 6.4 |
| 2 | 9 / 21 (SMT) | `--n-recv 128` | 6.45 | 7.32 | 0.46 | 7.0 |
| 4 | 5 / 9 (SMT) | `--n-recv 256` | 6.78, 6.78 | 5.72, 5.53 | 0.35 | 5.6 |
| 4 | 2 / 5 | `--n-recv 256` | 5.72 | 3.08 | 0.19 | 3.0 |

**The wire is at 12.1-12.2 GB/s** (`flood.c`'s 400 MB ping-pong over PSM2, two disjoint pairs, Open MPI's own choice
included; 97 % of the 12.5 GB/s of 100 Gb/s, and above session 19's 11.6-11.7 on the same environment).  **One rank a
host, the session 19 team, runs 16 hosts at 4.5-4.8 G points/s in PROBE and 4.0-4.1 G inserts/s in FILL: 0.28-0.30 and
0.25-0.26 a host, against 0.30-0.33 and 0.31 on 8 hosts** -- the scaling from 8 to 16 is linear within 10 %, and the
posted receives (`--n-recv 64`, 4 a peer for 15 peers) recover FILL's share (4.7, +15 %, robust: 4.55 at 128) without
touching PROBE.  The service thread is 90 % busy in PROBE, 72 % in FILL, the Router's verdict `service/network` in
every run: the limit of session 19, per byte the one MPI thread of the node hands to the kernel's SDMA path.  **Two
ranks a host, one bound to each socket, give each host two service threads and the engine 7.1-7.3 G points/s in PROBE
(+50 %) and 6.3-6.5 G inserts/s in FILL (+55 %) with the 8/22 team a rank and `--n-recv 128`** (three runs within
2 %): **0.45 a host, the same rate as one grvingt node alone at 10/21 (session 17's 0.46)**, the wire at 7.0 GB/s a
host (57 % of it), the service 85 % busy in PROBE and 70 % in FILL.  9/21 is 7.32 on one run, 10/20 6.1-6.9 and 7/23
6.6 (producer-bound: its senders wait 6 % of the phase): 8 or 9 producers a socket feed 21-22 dict threads there.  The one-a-core team (5/10 a rank) is
the FILL winner at 6.45-6.56 but loses a quarter of PROBE to a pool starvation the posted receives do not cure (below).
The SMT team on one rank (24/38) loses 5-10 % on both phases, as on gros (session 23).

## 1. Getting 16 nodes up

Nothing new in the recipe (session 19's `post19.sh` became `post24.sh`, one pass a node as root: host key forgotten,
`iommu.passthrough=1` in the cmdline, `Default domain type: Passthrough`, `rdma dev set opap94s0 name hfi1_0`, link
`ACTIVE` at 100 Gb/s, `cache_size` 4096, `numa_balancing` 0), and one thing to know: **the kadeploy3 client lost its
server connection after 855 s** (`Invalid request on kadeploy.nancy.grid5000.fr:25300 (Net::ReadTimeout)`, exit 1, an
empty log) **while the deployment itself had succeeded on all 16 nodes** -- every node answered as root on the new
kernel with 9 minutes of uptime.  Check the nodes, not the client's exit status.  The job's nodes were up 5 minutes
after `oarsub`, deployed 15 minutes after, benchmarking 18 minutes after.  The merge of `omp_reboot` (`0bff8a2`) into
`bench-grvingt` was clean and the release build of the demo took 20 s on the node.

## 2. The wire

`flood24.sh` (session 19's, `mpicc -O2 flood.c`, `mpiexec -n 2 --map-by ppr:1:node`, three runs a transport):

| pair | `--mca pml cm --mca mtl psm2` | Open MPI's choice |
| --- | --- | --- |
| grvingt-1 <-> grvingt-2 | 12.23, 12.10, 12.09 GB/s | 12.23, 11.75, 12.23 |
| grvingt-25 <-> grvingt-27 | 12.20, 11.05, 12.21 | 12.06, 11.94, 12.21 |

12.1-12.2 GB/s, 97 % of the line rate, on both pairs and with both selections (Open MPI picks PSM2 by itself once the
device is named `hfi1_0`).  Session 19 measured 11.6-11.7 on the same environment and the same program: the
difference is within what the two odd runs (11.05, 11.75) show, and the number to keep is **12.2 GB/s a direction**.

## 3. One rank a host: 16 hosts against 8

| P / R | knobs | FILL s | FILL G/s | PROBE s | PROBE G/s | wire a host, FILL / PROBE | service busy, FILL / PROBE | receives left unposted, FILL / PROBE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 / 21 | -- | 9.8 | 4.08 | 113.5 | 4.84 | 3.8 / 4.5 GB/s | 72 / 90 % | 7.9 M / 12.0 M |
| 10 / 21 | -- | 9.7 | 4.13 | 121.5 | 4.52 | 3.9 / 4.2 | | |
| 10 / 21 | -- | 10.0 | 4.00 | 122.7 | 4.48 | 3.7 / 4.2 | | |
| 10 / 21 | `--n-recv 64` | 8.5 | 4.71 | 114.1 | 4.82 | 4.4 / 4.5 | | |
| 10 / 21 | `--n-recv 128` | 8.8 | 4.55 | 122.2 | 4.50 | 4.3 / 4.2 | | |
| 24 / 38 | -- | 10.4 | 3.85 | 125.7 | 4.37 | 3.6 / 4.1 | | |

Session 19's 8 hosts did 2.40-2.64 G points/s in PROBE and 2.50 in FILL with the same team and the same collision
rate.  Twice the hosts give 1.7-2.0x on PROBE (the 8 % run-to-run scatter of the fabric is in both) and 1.6x on FILL --
1.9x with the posted receives at 4 a peer.  The Router's verdict is `service/network` in every phase of every run, the
service thread 72 % busy in FILL and 90 % in PROBE, the senders waiting 28-38 % of the phase and the receivers 22-50 %:
the node's one MPI thread is the ceiling, as in session 19, and the wire idles at 4.5 of its 12.2 GB/s.  The posted
receives matter to FILL only (every host sends every point off-node, and the default 32 slots are two a peer for
15 peers), by 12-15 %, and the SMT team on one rank loses: the second hardware thread of each core has nothing to
feed when the service does not keep up.

## 4. Two ranks a host, one a socket

| P / R a rank | knobs | FILL s | FILL G/s | PROBE s | PROBE G/s | wire a host in PROBE | service busy, FILL / PROBE | senders / receivers waited in PROBE | receives left unposted in PROBE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5 / 10 | -- | 6.2 | 6.45 | 156.5 | 3.51 | 3.4 GB/s | 64 / 29 % | 51 / 41 % | 188 M |
| 5 / 10 | -- | 6.2 | 6.45 | 155.6 | 3.53 | 3.4 | 64 / 29 | 51 / 41 | 438 M |
| 5 / 10 | `--n-recv 128` | 6.2 | 6.45 | 156.2 | 3.52 | 3.4 | 70 / 31 | 51 / 41 | 95 M |
| 5 / 10 | `--n-recv 128` | 6.5 | 6.16 | 153.7 | 3.58 | 3.4 | | | |
| 5 / 10 | `--n-recv 256` | 6.1 | 6.56 | 156.4 | 3.52 | 3.4 | 72 / 32 | 51 / 41 | 118 M |
| 6 / 9 | -- | 7.9 | 5.07 | 146.6 | 3.75 | 3.6 | 53 / 30 | 57 / 30 | 1071 M |
| 6 / 9 | `--n-recv 128` | 7.9 | 5.07 | 141.9 | 3.87 | 3.8 | | | |
| 10 / 20 | -- | 8.5 | 4.71 | 85.6 | 6.42 | 6.2 | 44 / 74 | 33 / 26 | 30 M |
| 10 / 20 | -- | 8.5 | 4.71 | 87.4 | 6.29 | 6.0 | 45 / 72 | 35 / 28 | 44 M |
| **10 / 20** | **`--n-recv 128`** | **6.3** | **6.35** | **89.5** | **6.14** | 6.0 | 70 / 71 | 36 / 29 | 19 M |
| 10 / 20 | `--n-recv 128` | 6.1 | 6.56 | 80.3 | 6.85 | 6.6 | 73 / 80 | 28 / 21 | 5 M |
| 10 / 20 | `--n-recv 128` | 6.0 | 6.67 | 86.9 | 6.33 | 6.2 | 73 / 73 | 34 / 27 | 14 M |
| 10 / 20 | `--n-recv 128 --credit 8` | 8.2 | 4.88 | 83.5 | 6.58 | 6.4 | 50 / 75 | 31 / 25 | 36 M |
| 8 / 22 | `--n-recv 128` | 6.2 | 6.45 | 77.8 | 7.07 | 6.8 | 71 / 84 | 9 / 26 | 0.2 M |
| 8 / 22 | `--n-recv 128` | 6.4 | 6.25 | 76.0 | 7.23 | 7.0 | 69 / 85 | 7 / 25 | 0.3 M |
| 8 / 22 | `--n-recv 128` | 6.3 | 6.35 | 77.0 | 7.14 | 7.0 | 70 / 85 | 8 / 25 | 0.3 M |
| 7 / 23 | `--n-recv 128` | 6.9 | 5.80 | 83.3 | 6.60 | 6.4 | 66 / 74 | 6 / 36 | 0 |
| 9 / 21 | `--n-recv 128` | 6.2 | 6.45 | 75.1 | 7.32 | 7.0 | 71 / 86 | 15 / 19 | 0.5 M |

Two things happen at once.  **With two service threads a host, FILL runs at 6.4-6.6 G inserts/s on the one-a-core
team** (+58 % over one rank; the service 64-72 % busy, the senders waiting 10-17 %) -- FILL was the phase session 19
found service-bound, and doubling the thread that binds it pays in full.  **PROBE on the same team collapses to 3.5**,
and the counters say why it is not the service: it is 29-32 % busy, the senders wait 51 % of the phase, the receivers
41 %, and the receiving side reports 95-438 M receives left unposted for want of a free block -- a starved pool, with
15-19 M sealed blocks parked for credit on the sending side.  More posted receives (128, 256 for 31 peers) cut the
unposted count by 2-4x and move PROBE not at all: the phase is pool-bound, not slot-bound.  **The SMT team a rank,
10/20, does not starve**: 6.1-6.4 G points/s in PROBE (+35 % over one rank, 0.38-0.40 a host, the wire at 6.0-6.2 GB/s
a host), the service 70-74 % busy, 19-44 M unposted.  Its FILL at the default receive count is 4.7, below one rank's
`--n-recv 64`; **with `--n-recv 128` FILL is 6.35-6.67 and PROBE 6.14-6.85** (three runs; the fabric's scatter is
11 % on PROBE here).  `--credit 8` on top gives FILL back to 4.88 (the senders wait 47 % of FILL for a free block: more
credit is more of the pool parked) and leaves PROBE alone.  **Shifting two producers to the dict side, 8/22 a rank,
is the best PROBE of the session and the tightest: 7.07, 7.23, 7.14** (77 s for 2^39), FILL 6.25-6.45, the service
85 % busy, the senders waiting 7-8 % of PROBE and the receivers 25 %, 0.3 M receives unposted; 9/21 gives 7.32 and
7/23 6.60 with its producers saturated (5.6 % wait).  PROBE follows the dict-thread count a host as long as the
producers keep up: 40 threads 6.1-6.9, 42 threads 7.3, 44 threads 7.1-7.2, 46 threads 6.6 (producer-bound).
**Hypothesis** (labelled): the 5/10 rank has 10 inboxes of 64 blocks and 10 receivers draining them; the 10/20 rank has
20 of each with the same pool (6383 blocks a rank); when 31 peers each park blocks for credit, the rank with fewer,
slower-drained inboxes latches into the exhausted-pool state session 8 saw across hosts, and the one with twice the
drain does not.  `--credit`, `--inbox` and the pool's size (`--block`, `--swc`) are the knobs to price; not measured
here beyond `--credit 8` (below).

## 5. Four ranks a host, two a socket

`--map-by ppr:2:socket:pe=8 --bind-to core`, 10G a rank, 63 peers so `--n-recv 256`; each Router sees 16 CPUs (8
cores), keeps the service alone on one and has 7 cores for the workers: 2/5 one a core, 5/9 two a core.

| P / R a rank | FILL s | FILL G/s | PROBE s | PROBE G/s | wire a host in PROBE | service busy, FILL / PROBE | senders / receivers waited in PROBE | receives left unposted in PROBE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5 / 9 | 5.9 | 6.78 | 96.1 | 5.72 | 5.6 GB/s | 48 / 28 % | 41 / 29 % | 12 M |
| 5 / 9 | 5.9 | 6.78 | 99.4 | 5.53 | 5.6 | 47 / 27 | 43 / 31 | 19 M |
| 2 / 5 | 7.0 | 5.72 | 178.6 | 3.08 | 3.0 | 28 / 10 | 47 / 47 | 90 M |

Four service threads a host make **FILL the fastest of the session, 6.78 G inserts/s** (+65 % over one rank), with
zero receives left unposted in that phase and the services 48 % busy.  **PROBE is 5.5-5.7**, below the two-rank 10/20
team's 6.1-6.4 with fewer dict threads a host (36 against 40) and the four services at 28 % busy while the senders wait
41-43 % and the receivers 29-31 %: the same pool-bound state as the two-rank 5/10 team, milder.  The one-a-core team
(2/5 a rank, 20 dict threads a host) is the worst PROBE of the session at 3.08 with the services 10 % busy and 90 M
receives unposted.  The pattern over the three rank counts is one pattern: **PROBE follows the dict-thread count a
host (20-21 threads: 3.1-4.8; 36: 5.5-5.7; 40: 6.1-6.4) once each host has at least two service threads, and a
rank with few receivers under many peers starves its pool**; FILL follows the service-thread count a host (1: 4.0-4.7,
2: 6.3-6.6, 4: 6.8) until the wire or the pool takes over.

## What to change

1. **The 16-host grvingt reference is two ranks a host, one bound to each socket** (`--map-by ppr:1:socket --bind-to
   socket`, `-x HFI_NO_CPUAFFINITY=1`), **8 producers + 22 dict threads a rank on the socket's 15 free cores,
   `--ram 20G` a rank, `--n 39`, `--prefetch 8 --n-recv 128`: FILL 6.3-6.5 G inserts/s, PROBE 7.1-7.3 G points/s**,
   0.45 a host = one node's own rate, 57 % of the wire.  One rank a host with the session 19 team is 4.0-4.1 / 4.5-4.8,
   and its posted receives want `--n-recv 64` (4 a peer) for FILL.
2. **Ranks a host is the lever wherever the service thread binds**, on Omni-Path as on 25 GbE (session 22): a second
   service thread a host is worth +50 % in PROBE and +55 % in FILL here, for no code, and brings the per-host rate
   back to the single node's.  The Router's one-service-thread
   design is what makes the socket the unit; the alternative is a second service thread inside the Router.
3. **The exhausted-pool latch is back**, in the 5/10-a-rank configuration under 31 peers: PROBE at 3.5 with everyone
   waiting and the service idle.  It is the class of failure session 8 fixed across 4 hosts and session 22 relieved
   with `--n-recv`; here the receive slots do not cure it.  Price `--credit`, `--inbox` and the pool before trusting
   any small team over many peers.
4. **Posted receives scale with the peers**: 4 a peer (`--n-recv 64` for 15, 128 for 31) is worth 12-35 % of FILL and
   nothing on PROBE.  The default of 32 should become `4 * (n_nodes - 1)`, floored at 32.
5. **The kadeploy client's exit status is not the deployment's**: it timed out on the server after 855 s with all 16
   nodes deployed.  Check the nodes.

## Reproducing

```bash
# frontend: reserve, deploy (the client may time out; check the nodes), per-node setup, wire check
oarsub -t deploy -p grvingt --project cryptanalyse -l nodes=16,walltime=3:00 "sleep 10800"
bash build/session24/deploy24.sh        # waits for the job, writes nodes.txt, kadeploy3 -a ~/debian13-big-iommupt.yaml -f nodes.txt -k
bash build/session24/post24.sh          # per node as root: cmdline, Passthrough, OPA rename, link, hfi1 cache_size 4096, numa_balancing 0
bash build/session24/flood24.sh         # 12.1-12.2 GB/s over PSM2 between the first two nodes
# grvingt-1, login shell: build, one run, the batches
ssh grvingt-1 bash -l build/session24/setup24.sh
P=10 R=21 bash -l build/session24/run24.sh                                   # one rank a host, --n 39 --ram 40G
PPN=2 P=8 R=22 EXTRA="--prefetch 8 --n-recv 128" bash -l build/session24/run24.sh    # the reference, a rank a socket, 20G a rank
PPN=2 RAM=20G CONFIGS="a|10|20|--n-recv 128;b|5|10|" bash -l build/session24/batch24.sh
python3 build/session24/summ24.py ref1 recv64 p2_10_20_r128 p2_8_22_r128 ...
```

# Session 26 -- grvingt, PCS: what makes a search finish, not what makes a round fast

2026-09-16/17, **grvingt-12** (2 x Xeon Gold 6130, 32 cores / 64 PU, 187 GB), `oarsub -p grvingt --project
cryptanalyse -l host=1,walltime=6:00:00` (job 6929764), standard `debian13` environment, kernel 6.12.107,
`module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0 gcc-toolchain`, `kernel.numa_balancing = 0`,
transparent huge pages `[always]`.  Commit `0f8a000` on `bench-grvingt` (= `0bff8a2` of `omp_reboot`, the placer
keeping the service thread alone on its core), release build, the one session 25 used.  One rank, `--engine pcs`,
`double_speck64_demo`, seed 1 everywhere, 61 workers + 1 service in every run.  Scripts and logs in
`build/session26/`; `summ26.py` is the table below.

**The question** (the user): get the *result* as fast as possible, i.e. minimise `#rounds * t_round`, over the
dictionary size w (up to 160 GB) and the team split -- all the RAM need not be used, and a harder function might
pay by cutting the memory traffic.

## The figure of merit

A round is worth the **distinct collisions it finds**, `E_i` (the engine prints it, `#distinct coll (this i)`,
from the walkers' merged HyperLogLog over the problem-level pairs).  For a claw problem whose f and g are random
functions of an n-bit domain, the pairs a run can ever record -- f-f, g-g and f-g -- number `U = 2^(n+1)`, a
round's located set is a uniform sample of those active under its version of the mixing, and the golden pair is
one of them, so

    P(found in a round) = E_i / 2^(n+1)     #rounds = 2^(n+1) / E_i     total = 2^(n+1) * (t_round / E_i)

**`t_round / E_i`, the wall-clock per distinct collision, is therefore the whole objective**, and it is what the
`ns/coll` column below is.  w drops out of it: a dictionary twice the size makes a round twice as long and worth
twice as many collisions.  Two rounds of a configuration measure it; nothing has to run to completion.

It splits into two factors that move in opposite directions, which is the whole story of this session:

    ns/coll  =  ev/coll / rate       ev/coll = k * beta / (theta * C)      C = E_i / w

`ev/coll`, the evaluations per distinct collision, is **algorithmic and machine-independent**; `rate` is the
node's evaluations per second.  `C` depends on **alpha alone** (1.04 at alpha 2.5 whether n is 30, 36 or 42 and
whether w is 12.5 M or 20 G slots), `1/theta = alpha^-1 sqrt(2^n / w)` depends on the dictionary's share of the
domain, and `k` (evaluations over pure walking) is the resolver's overhead.

## Summary

| | |
| --- | --- |
| Best on the node, 160 GB, n=42 | `--ram 160G --producers-per-node 39 --dicts-per-node 22`, 52.6 ns/coll |
| Difficulty: optimum | alpha 1.8-2.5 (`1/theta` 4.7-6.5), flat to 1 %; the shipped default is in it |
| Difficulty: alpha 5 (twice the dictionary's worth) | 2.0x worse -- `C` collapses from 1.04 to 0.38 |
| Difficulty: alpha 0.35 | 3.2x worse -- 262 evaluations a collision against 51 |
| RAM: 160 GB vs 40 / 10 / 2.5 GB, n=42 | 1.00 : 1.36 : 2.12 : 3.87.  **Use every byte** |
| RAM: cost of a big dictionary per evaluation | none: 30.4 M f/s a walker at 160 GB, 30.2 at 4 GB |
| beta | 8 is the optimum; 4 and 16 cost 8 %, 2 costs 40 %, 32 costs 26 % |
| Split | `R` dict threads follow the DP rate: `R ~ 220 * theta / (1 + 3.6 * theta)`, `P = 61 - R` |
| Round count | `2^(n+1) / E_i` confirmed over 12 full searches at n=30 (1.29x the model, SE 0.17) |

## 1. The difficulty against the split, at a fixed dictionary

`--n 36 --ram 4G` (w = 500.0 M slots = 2^28.90), `--alpha` swept, 2 rounds each, **2-3 splits per alpha and not
an exhaustive search of them** -- three around the bottom of the U, two at its ends, the best of each below.  `1/theta` is what alpha buys at this w; the round is `beta * w` distinguished points either way.

| alpha | 1/theta | P/R | round | ev/coll | node eval/s | DP/s | dist/w | **ns/coll** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5.00 | 2.34 | 31/30 | 18.1 s | 57.8 | 0.61 G | 221.8 M | 0.38 | 95.1 |
| 3.50 | 3.35 | 31/30 | 19.1 s | 51.0 | 0.90 G | 209.7 M | 0.67 | 56.9 |
| 2.50 | 4.69 | 36/25 | 24.5 s | 51.0 | 1.08 G | 163.7 M | 1.04 | 47.0 |
| **1.80** | **6.51** | **41/20** | 33.8 s | 56.5 | 1.21 G | 118.5 M | 1.45 | **46.5** |
| 1.30 | 9.02 | 46/15 | 47.3 s | 68.1 | 1.33 G | 84.7 M | 1.85 | 51.0 |
| 0.90 | 13.03 | 50/11 | 68.5 s | 91.0 | 1.47 G | 58.4 M | 2.21 | 61.9 |
| 0.60 | 19.54 | 55/6 | 102.0 s | 135.5 | 1.62 G | 39.2 M | 2.43 | 83.8 |
| 0.35 | 33.50 | 57/4 | 178.4 s | 262.1 | 1.74 G | 22.4 M | 2.37 | 150.3 |

**A harder function does buy throughput, and it is not enough.**  The node's rate climbs monotonically with the
difficulty, 0.61 -> 1.74 G evaluations/s, exactly as hoped: fewer distinguished points per evaluation, less
dictionary traffic, dict threads turned into walkers.  But the walking itself grows faster: `ev/coll` bottoms out
at 51 around alpha 2.5-3.5 and is 5x that at alpha 0.35.  The product is flat over alpha 1.8-2.5 and rises on
both sides.

**The easy side is a cliff, the hard side a slope.**  At alpha 5 the round is 26 % shorter than the default's
but finds 0.38 *w* distinct collisions instead of 1.04 *w*: the trails are too short for the dictionary to be
worth filling (`1/theta` = 2.34, `beta * w` points arriving at w slots off trails of 2.3 steps).  That is the
regime a bigger dictionary at a fixed alpha would walk into, and it costs 2x.

## 2. The RAM question: one problem, four dictionaries

`--n 42`, alpha at its 2.5 default, each size at the best of 2-3 splits.  This is the user's question in its
own terms: 160 GB is what the node holds, and the alternative is a smaller dictionary and a harder function.

(160 GB and 40 GB are one round each, so their times carry the shard's first touch; 10 GB and 2.5 GB are the
second of two.  The first touch is parallel over the dict threads, seconds against a round of minutes.)

| RAM | w | 1/theta | P/R | round | ev/coll | node eval/s | **ns/coll** | vs 160 GB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **160 G** | 20.0 G | 5.93 | 39/22 | 1125.3 s | 62.3 | 1.19 G | **52.6** | 1.00 |
| 40 G | 5.0 G | 11.86 | 46/15 | 407.2 s | 115.6 | 1.62 G | 71.4 | 1.36 |
| 10 G | 1.25 G | 23.73 | 52/9 | 161.7 s | 226.6 | 2.03 G | 111.5 | 2.12 |
| 2.5 G | 312.5 M | 47.45 | 56/5 | 75.1 s | 461.7 | 2.27 G | 203.6 | 3.87 |

**Every byte earns its place.**  Distinct collisions per round stay at 1.07-1.18 *w* across the four, so the
round count falls exactly as w rises, while `ev/coll` doubles with every 4x cut (62.3, 115.6, 226.6, 461.7 --
the `1/sqrt(w)` of the theory, to 3 %).  The rate does climb, 1.19 -> 2.27 G/s, and it pays for a third of the
loss at the first step and less after: net 1.36x, 2.12x, 3.87x.  There is no size at which stopping short wins.

The other reading of the same table: at 160 GB neither side of the node is starved.  Its walkers wait 5.3 % of
the round and its dict threads 5.4 %, both under `Router_Bottleneck`'s 5 % threshold -- which is why the verdict
line falls through to `service/network`, a meaningless label at np 1 -- and the neighbouring split, 36/25, leaves
the dict threads idle 18 % and is slower (53.3 ns/coll).  Near-balance is an inference from those two numbers,
not something the Router says.  Below 160 GB the split is what recovers most of the difference: 40 GB at the
predicted 50/11 was 88.7 ns/coll, at 46/15 it is 71.4.

## 3. A big dictionary costs nothing per evaluation

The control: the same difficulty as the 4 GB baseline with 16x the dictionary, `--n 40 --ram 64G`, split 36/25.

| | round | ev/coll | node eval/s | f/s a walker | probe/s a dict thread | ns/coll |
| --- | --- | --- | --- | --- | --- | --- |
| 4 GB, n=36 | 24.4 s | 51.1 | 1.09 G | 30.2 M | 6.6 M | 46.8 |
| 64 GB, n=40 | 403.4 s | 50.9 | 1.05 G | 30.0 M | 6.3 M | 48.4 |
| 160 GB, n=42 (1/theta 5.93) | 1125.3 s | 62.3 | 1.19 G | 30.4 M | 6.5 M | 52.6 |

16x the dictionary costs 3.7 % of the rate, 40x costs nothing measurable: a dict thread does 6.3-6.6 M probes/s
and a walker 30 M evaluations/s whatever the shard's size.  Transparent huge pages are `[always]` on this image
and `PcsDict`'s `vector<u64>` gets 2 MB pages at every size, so the TLB reach does not change with w.  **Every
cheap 4 GB measurement therefore transfers to the full-RAM scale**, which is what makes section 1 usable.

## 4. beta, the round's length

`--n 36 --ram 4G`, alpha 2.5, split 36/25 throughout (tuned for beta 8, not re-tuned per beta).

| beta | round | ev/coll | dist/w | ns/coll |
| --- | --- | --- | --- | --- |
| 2 | 5.6 s | 67.0 | 0.17 | 65.8 |
| 4 | 11.4 s | 54.8 | 0.45 | 50.6 |
| **8** | 24.5 s | 51.0 | 1.04 | **47.0** |
| 16 | 52.5 s | 54.1 | 2.07 | 50.6 |
| 32 | 110.3 s | 62.7 | 3.73 | 59.0 |

The default is the optimum and the curve is shallow around it.  Short rounds waste the dictionary (at beta 2 it
is 18 % full when the round ends); long ones waste the arrivals, since a slot keeps one trail and `PcsDict`'s
most-recent-wins policy throws the rest away.

## 5. The round count itself

The other half of `#rounds * t_round` is not visible in any single round, so: 12 full searches, `--n 30
--ram 100M` (w = 12.5 M slots, `1/theta` 3.71, `E_i` = 1.00 *w* every round), 36/25, no `--nrounds` cap, seeds
1-12.  Rounds to the golden pair: 10, 57, 81, 168, 184, 221, 227, 313, 356, 361, 471, 635 -- **mean 257** against
`2^31 / E_i` = 172 predicted (a geometric, so the 12-seed standard error is 74: +1.2 sigma).  The sharper
statistic is the coverage at the find, whose mean is `U/2 = 2^n` with a standard error of `U/(2*sqrt(3*12))`:
measured **1.387 G distinct pairs against 1.074 G predicted, a ratio of 1.29 +- 0.17**.

So `#rounds = 2^(n+1) / E_i` is the right form and **the model holds to within 30 %** -- but both statistics
land high (+1.2 and +1.75 sigma), which is more likely a small systematic than noise: a planted pair slightly
harder than a random collision, or located collisions favouring long trails.  Read every absolute below with
that margin.  One line so that the exponent is not re-litigated: `U` is `2^(n+1)` rather than `2^n` because f
and g are **random functions** (Speck under a key, truncated), not permutations; a problem built from two
permutations would record `2^(n-1)` pairs, all of them claws.  The n=30 run is what settled it.  It also says a
search is *half the problem's collisions* long: every PCS run ends deep in the regime where rounds re-find pairs
earlier rounds already had, and the geometric above is what accounts for it.

**Absolute wall times on this node** (`total = 2^(n+1) * ns/coll`): n=36 with 4 GB, 1.8 h.  n=42 with 160 GB at
39/22, **5.3 days** -- 5 to 7 with the margin above.  n=48 against that same 160 GB would be 2^9 times it -- the domain grows by 2^6 and
`1/theta` with it by 2^3 -- i.e. 7.4 years on one node, which is what distributing w is for.

## 6. The split, as a rule

A dict thread retires 5.3-7.5 M probes/s and the optimum leaves it a little headroom, so the dict thread count
follows the node's distinguished-point rate, not the difficulty as such:

    R ~ 220 * theta / (1 + 3.6 * theta)     P = (workers) - R

which gives 26, 22, 17, 14, 13 at `1/theta` = 4.69, 6.51, 9.02, 11.86, 13.03 against the 25, 20, 15, 15, 11
measured -- close enough to seed a two-point sweep.  The operational version is simpler: **add dict threads until
the senders stop waiting**, which the round report's `LIMITED BY` line names outright.  Getting it wrong is
expensive: 40 GB at n=42 costs 24 % between 46/15 and 50/11, and 71 % between 46/15 and 53/8.

## Where this sits against earlier sessions

Session 25's baseline reproduces exactly: 4 GB, n=36, 36/25 in 24.4 s here against 24.5 s there, 164.6 M DP/s,
1.04 *w*.  Its "the real difficulty is dictionary-bound, the easy one walker-bound" stands, and this session
says what to do about it: nothing -- the point where the dictionary starts to bind is also where the algorithm
is cheapest, and the two optima coincide to within 1 %.  Session 25's own entry was never written into this file
(its reservation ended first); its numbers survive in `build/session25/` and in the project memory.

Nothing here contradicts the direct engine's sessions.  One prediction this session does *not* measure: at
alpha 1.8 the node ships 118 M DP/s instead of 164 M for the same work, **28 % less wire traffic per evaluation**,
which is free on one node and should not be on the 8-16 host runs where the wire or the service thread was the
limit (sessions 20, 24).  That wants a multi-node A/B before it is believed.

## What to change

- **Nothing in the code.**  `--alpha 2.5` and `--beta 8` are both at their optimum, and `--ram` should be
  everything the node has.
- On one grvingt node, the PCS reference is `--ram 160G --producers-per-node 39 --dicts-per-node 22` at n=42,
  or `R` from the rule of section 6 at another size.
- `--alpha 1.8` is worth trying on a multi-node run, where its 28 % less wire traffic per evaluation is not free.

## Reproducing

```bash
oarsub -p grvingt --project cryptanalyse -l host=1,walltime=6:00:00 'sleep infinity'
# on the node
sudo-g5k sysctl -w kernel.numa_balancing=0
cd ~/mitm-grvingt/build/session26
bash sweepA.sh          # difficulty x split, 4 GB, n=36          ~40 min
bash rest.sh            # the RAM question, the control, the round count, beta   ~2 h 30
python3 summ26.py a_*.log        # one line per run, sorted by ns/coll
```

One run alone is `P=39 R=22 N=42 RAM=160G ROUNDS=1 bash run26.sh` (`EXTRA="--alpha 1.8 --beta 8"` for the rest).
The figure of merit is `secs / (dist/w * w)`; `summ26.py` does it, and `w` is the banner's "slots in all".
