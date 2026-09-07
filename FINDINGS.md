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
