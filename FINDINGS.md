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

# Session 14 -- grdix on HDR200: the default block costs 38 %, and one thread is the reason

2026-09-12, **grdix-8** and **grdix-9** (2 x EPYC 9754, Zen 4c, 256 cores / 512 PU, 2 NUMA nodes,
32 L3 of 16 MB, 1 TB), `oarsub -I -p grdix --project cryptanalyse --resource host=2,walltime=2`
(job 6922626).  Commit `23e713b` on `bench-grdix`, release build (`-O3 -DNDEBUG -march=native
-ggdb`).  The image is now **debian13** and carries no MPI: `module load openmpi/4.1.6` gives a
Guix-built Open MPI **4.1.6** under `/gnu/store`, and hwloc is a module too (no `/usr/include/
hwloc.h`), so CMake needs `HWLOC_ROOT` pointed at the Guix prefix.  Interconnect: **Infiniband
HDR200** (`ibp66s0`, 200 Gb/s, ACTIVE) alongside 25 GbE on `br0` and two RoCE ports; `--mca pml
ucx` with `UCX_NET_DEVICES=ibp66s0:1`.

**These numbers do not compare to sessions 4, 6 or 11.**  `45f59ff` made `router_bench` a *pure*
router benchmark: the sender's PRNG is gone (a counter spread by a Fibonacci multiply picks the
destination, no modulo) and the receiver's `murmur128` is gone (an add and a XOR).  The ~12 ns of
PRNG that used to sit on every push is what those sessions' 2-3 G points/s were mostly measuring.
The accounting check (`pushed != popped`) stayed silent throughout: nothing was lost in any run
below.

**`xor` is `0000000000000000` in every round, and it is a dead check in this bench.**  Not a bug
and not sender-count parity: a sender pushes `(n, n)` for `n = 0, 1, 2, ...` and stops on a
multiple of 1024, a receiver folds `a + b = 2n`, and the XOR of `2n` over any 1024-aligned run of
`n` is exactly zero -- the low ten bits cancel pairwise and every higher bit appears 1024 times.
So every sender contributes 0 whatever it does, and the fold cannot catch a lost or duplicated
block.  Only the `pushed != popped` line can.

## Summary

| | one node (grdix-8) | two nodes (grdix-8 + grdix-9) |
| --- | --- | --- |
| default block (4096 points) | 8.1 G points/s | 1.5 G points/s per node |
| block 16384 | **13.0 G points/s** | **2.5 G points/s per node** |
| what limits it there | memory: 416 of ~442-561 GB/s | fabric: 19.96 of 23.54 GB/s each way |
| what limits the default | a per-block cost, ~2.0 M blocks/s (section 6) | the same, then the fabric |

## 1. The default block size costs 38 % on one node

`--senders 150 --receivers 105`, 2 rounds of 10 s, round 1 reported.  The pool holds a **fixed
number of blocks** (16808 for this team), so `block_points` scales its memory linearly.

| `--block` (points) | message | pool | routed | ns/point | blocks/s |
| --- | --- | --- | --- | --- | --- |
| **4096 (the default)** | 65600 B | 1.10 GB | 8.1 G | 18.6 | 2.0 M |
| 8192 | 131136 B | 2.20 GB | 12.7 G | 11.9 | 1.6 M |
| 16384 | 262208 B | 4.41 GB | **13.0 G** | 11.7 | 791 K |
| 32768 | 524352 B | 8.81 GB | 12.9 G | 11.6 | 395 K |
| 65536 | 1048640 B | 17.63 GB | 13.0 G | 11.6 | 198 K |

The knee is between 4096 and 8192 and the curve is flat from 8192 up.  Nothing else was changed.

## 2. The ceiling is a per-block cost, at ~2.0 M blocks/s

The number that does not move is **blocks per second**.  Across every team size from 128 workers
up, every sender/receiver split, every `--sweep` and every `--inbox`, the default block sits at
**2.0 M blocks/s** and the point rate is just that times `block_points`:

| config (block 4096, 150 + 105 unless stated) | routed | blocks/s |
| --- | --- | --- |
| default (`--sweep 256 --inbox 64`) | 8.1 G | 2.0 M |
| `--sweep 1024` | 8.3 G | 2.0 M |
| `--sweep 4096` (916 turns/s, 2183 blocks a turn) | 8.2 G | 2.0 M |
| `--inbox 256` (36968 blocks, 2.43 GB) | 8.6 G | 2.1 M |
| `--inbox 768` (90728 blocks, 5.95 GB) | 8.9 G | 2.2 M |
| 128 + 127 / 192 + 63 / 224 + 31 | 8.0 / 8.3 / 9.6 G | 1.95 / 2.03 / 2.3 M |

So it is neither the sweep cap (raising it 16x moves nothing, and at `--sweep 4096` the service
handles 2183 blocks a turn without going faster) nor the buffer depth (`--inbox 768` gives 5.4x
the pool -- **more** points in flight than `--block 16384` has -- for 10 %).  It is a fixed cost
per block on something serialized, and **the service thread is the only single-threaded thing in
the path**: every sealed block, local ones included, passes through it.

`perf record -C 0` (the placer pins the service to CPU 0) during a block-4096 round, 12 s,
`--sort srcline`:

| share of the service's cycles | source |
| --- | --- |
| 24.7 % + 21.0 % | `service.hpp:43-44` -- the `n_valid` scan inside `Router_node::complete()` |
| 16.0 % | `atomic_base.h:501` -- the `load_acquire` that scan is made of |
| 6.5 % | `service.hpp:249` |
| 6.5 % | `ring.hpp:48` |

`complete()` asks "are all L lines of this sealed block written?" by reading L `n_valid` bytes --
one cache line, since `nv_stride` rounds L up to 64 -- and that line is being written by the
senders at that very moment, so each read is a coherence miss, and an incomplete block is scanned
again on the next turn.  **The cost is the miss, not the length**: `--swc 1024` halves L to 4 and
changes nothing at all (8.0 G, the same 2.0 M blocks/s).  (`--swc 4096`, i.e. L = 1, is not the
other end of that line -- it drops to 4.1 G and 1.0 M blocks/s, because a sender must then fill a
64 KB private line before publishing anything.  Different effect, not a scan that got cheaper.)

The core is **not idle**: perf counted 26.26 G cycles in 12 s on CPU 0,
i.e. 2.19 GHz on a 2.25 GHz part, and ~46 % of them are in that scan -- about 500 cycles a block
on the scan alone.  (The round line's `service 0 % busy` is
`(TURNS - IDLE_TURNS)/TURNS`, which counts MPI work only.  At np 1 there is none, so it reads 0 %
while the thread is saturated.  That metric is misleading and should be renamed or fixed.)

The same profile at `--block 16384` is the control: `complete()` is **gone from the top of the
profile** (nothing above 5.3 %, and the list is `PMPI_Testsome` and
`ompi_request_default_test_some` -- the service polling an empty network), because there are 4x
fewer blocks to scan.

**That profile does not prove the service is the limiter, and section 6 shows it is not.**  The
service thread spins, so it burns 100 % of its cycles whatever happens; the share of them spent
in `complete()` says where its time goes, not that it is on the critical path.  Taking the scan
off it entirely buys 2.5 %.  What survives from this section is the *measurement*, which is solid
and is the thing to explain: a fixed cost per block, ~2.0 M blocks/s, immovable by `--sweep`,
`--inbox` and the split alike, so that the point rate is that rate times `block_points`.

## 3. Above the knee, the node is at its memory system

A point is 16 bytes, written into a block with non-temporal stores (so it lands in DRAM, never in
the writer's cache) and read back by a receiver.  13.0 G points/s is therefore **208 GB/s written
+ 208 GB/s read = 416 GB/s**.  An OpenMP STREAM on the same node (`~/membw.c`, 256 threads,
8 GB arrays, parallel first touch):

```
triad   :  498 GB/s   (counts the read-for-ownership: 4 streams)
nt-write:  442 GB/s   (_mm512_stream_pd, no RFO)
read    :  561 GB/s
```

416 GB/s of a half-read/half-write mix against 442 pure-write and 561 pure-read is the wall, not
a coincidence.  So the router has **two** ceilings on one node and the block size chooses which
one it meets: the service thread below 8192 points, the memory system above.

That also explains why the sender/receiver split looked irrelevant at the default block and is not
(np 1; the block-4096 row is 3 s rounds, the block-16384 row 6 s):

| split | 128 + 127 | 150 + 105 | 192 + 63 | 224 + 31 |
| --- | --- | --- | --- | --- |
| block 4096 | 8.0 G | 8.1 G | 8.1-8.3 G | -- |
| block 16384 | 12.2 G | **12.9 G** | 10.0 G | 9.6 G |

At the default block every split is pinned to the same 2.0 M blocks/s, so the split cannot show.
Once the service is out of the way, receivers earn their cores: `192 + 63` loses 22 % against
`150 + 105`.  Worker scaling at the default block, for the record: 8+7 1.7 G, 16+15 2.0 G,
32+31 3.2 G, 64+63 5.9 G, 96+95 7.2 G -- the last two are already within 12 % of the 2.0 M cap.

## 4. Two nodes: the router reaches 85 % of the wire, in both directions at once

`-H grdix-8,grdix-9 --npernode 1`, 150 + 105, 10 s rounds (session 9's minimum: the drain adds
1.4-3.5 s).  Destinations are uniform, so half of every node's points cross the link.  `net` is
bytes **sent** per node, one direction.

| `--block` | routed/node | net per node | msgs/s | service "busy" |
| --- | --- | --- | --- | --- |
| 4096 (default) | 1.5 G | 11.7-12.3 GB/s | 178 K | 39 % |
| 8192 | 2.2 G | 17.44 GB/s | 133 K | 18 % |
| 16384 | **2.5 G** | **19.96 GB/s** | 76 K | 6 % |
| 32768 | 2.5 G | 19.88 GB/s | 38 K | 2 % |
| 65536 | 2.5 G | 19.92 GB/s | 19 K | 1 % |

`192 + 63` at block 16384 gives the same 2.5 G and 19.82 GB/s, so the split stops mattering once
the fabric is the limit.  The reference is `ucx_perftest -t tag_bw -s 262144` between the two
nodes: **23.54 GB/s** one way, and **23.54 GB/s in each direction with both running at once** --
the link is honestly full duplex and the HCA's PCIe is not in the way.  The router sends 19.96
and receives 19.96 simultaneously: **85 % of the raw rate, both ways at once.**

The control is `--local-only` at np 2 (every point stays on its node): 12.8 G/node, i.e. exactly
the single-node rate of §1.  So running two ranks costs nothing in itself; the 2.5 G is the remote
half meeting a narrow pipe.

## 5. What this means

Compared consistently -- bytes one way against bytes one way -- a node writes **208 GB/s** into
blocks locally and can send **23.5 GB/s** to its peer: the fabric is about **9x narrower than
memory** (equivalently 13.0 G points/s local against 1.25 G points/s remote per node, ~10x).  So
a second node does not add per-node throughput, it divides it: 13.0 G alone, 2.5 G each in a
pair, because half of every node's points now take the narrow path.  In aggregate that is
**13.0 G points/s on one node against 5.0 G on two** -- at a uniform fan-out the second node makes
the whole thing 2.6x slower.  That is a property of the placement, not of the router: the router
is at 85 % of the wire, and the wire is what it is.

For the engines this is the number that matters: **a second node buys 2x the dictionary at about
1/5 the per-node rate.**  Nothing in the router can change that, and neither can smarter sharding
-- a point's shard is its hash, so the fan-out is uniform by construction.  What the measurement
does say is that the trade is worth making only when the dictionary has to be that big.

## 6. The experiment: the service thread is not the ceiling

Commit `762c660` moves the wait off the service and onto the sealer, in six lines: the sealer has
just installed the fresh block and blocks nobody, so it can spin on `complete()` itself and seal
only a whole block.  `handle_block` then dispatches unconditionally and the `pending` list is
dead.  The scan still happens -- but on as many sender cores as there are destinations being
sealed, in parallel, instead of on one thread in series.  If the service's scan were the ceiling
this had to lift it.

| `--block 4096`, 150 + 105, interleaved | routed | ns/point | blocks/s |
| --- | --- | --- | --- |
| baseline `63be9c4` (rep 1 / rep 2) | 8.1 / 8.1 G | 18.6 / 18.6 | 2.0 M |
| sealer waits `762c660` (rep 1 / rep 2) | 8.3 / 8.3 G | 18.0 / 18.1 | 2.0 M |

**+2.5 %, and blocks/s does not move at all.**  `router_test --test all` passes at np 1 and np 2
on the node.  So the service thread was never the constraint, and the 46 % of its cycles in
`complete()` was a spinning thread's idle time wearing a costume.

A whole-machine profile on the experiment build says the same thing from the other side.  At
`--block 4096` the cycles are spread over real work -- `Router_Push`'s body (workers.hpp:160-165,
~17 %), the non-temporal copy (common.hpp:90, 7.2 %), `Router_Pop` (workers.hpp:243-245, ~12 %) --
with no hotspot and, notably, **no receiver idling**.  At `--block 16384` the top line is
`workers.hpp:238` at **22 %**: that is `Router_Pop` finding no block, i.e. receivers starved
because the senders cannot feed them any faster.  So the two block sizes are limited at different
ends -- senders at 4096, the memory system at 16384 -- and what the big block relieves is
something on the *sender's* side of a block transition, not the service's.

## 7. It is not work, it is stall: three more mechanisms ruled out, and the counter that settles it

Sections 2 and 6 left "a fixed cost per block" unlocated.  Three interventions, each built on a branch
and A/B'd interleaved against `63be9c4` on the same node, at `--block 4096`, 150 + 105:

| intervention | what it removes | routed | blocks/s |
| --- | --- | --- | --- |
| baseline `63be9c4` | -- | 8.1-8.2 G | 2.0 M |
| `762c660` the sealer waits | the service's `complete()` scan | 8.3 G | 2.0 M |
| `c256db5` one sealed stack per group | the node-wide `sealed_top` CAS | 8.2-8.3 G | 2.0 M |
| `--dests 210` / `420` | destination-word transitions, 2x and 4x rarer | 7.9 / 8.0 G | 1.9 / 2.0 M |

and two more from the flags: with the senders held at 150, blocks/s *rises* with more receivers
(24 -> 1.8 M, 48 -> 1.9 M, 105 -> 2.0 M), which is the opposite of contention on the free ring's
ticket word; and `--swc` 128 / 256 / 512 gives 7.8 / 8.0 / 8.1 G, so the per-line cost is small and
additive, not a cliff.  **Every node-wide singleton on the sender's per-block path is now ruled out
by an intervention rather than by a profile.**

`perf stat -a` over 5 s of a steady round says why they all failed, and it should have been the
first measurement of the session:

| | `--block 4096` | `--block 16384` |
| --- | --- | --- |
| cycles | 4.013e12 | 4.020e12 |
| instructions | 2.074e12 | 3.238e12 |
| **instructions per point** | **51.2** | **49.4** |
| **IPC** | **0.52** | **0.81** |
| cache misses per point | 0.070 | 0.051 |

The cycles are equal because every thread spins.  The **instructions per point are equal too** --
within 3.6 % -- so the small block does not make the router execute more code.  It makes the same
code stall: IPC falls by 36 % and the machine takes 37 % more cache misses per point.  No lock, no
CAS and no extra branch can produce that signature, which is exactly why removing three of them
changed nothing.

It is not bandwidth either.  At `--block 4096` the node moves 131 GB/s written + 131 read = 262 GB/s,
against the 416 GB/s it sustains at 16384 and the 442-561 GB/s of section 3.  So the small block
leaves both the memory system *and* the cores idle at once: what it costs is **memory-level
parallelism**.  Converting the gap, the small block costs about 104 extra cache misses per block
transition, and none of them are in code that a profile attributes to a hot line.

**What to measure next**, and this is a measurement, not another guess: the miss distribution by role
and address.  `perf record -e cache-misses -C <sender cpus>` against `-C <receiver cpus>` says which
side stalls, and `perf mem record` / `perf c2c` says whether the misses are cold DRAM reads, the
non-temporal stores' own write-allocate behaviour, or a dependent chain in the install path.  The
one structural suspect that survives -- untested because the pool is sized by a formula rather than
a flag -- is the pool's *depth in points*: `n_blocks` is fixed at 16808 by the team's shape, so the
small block holds 68.8 M points in flight against 275 M, and section 2's `--inbox` test did **not**
rule that out (raising `inbox_blocks` grows the pool and the receivers' commitment by the same
`R * (inbox + 1)`, leaving the freely circulating count near 5200 either way).  Adding a slack knob
is a two-line change and the right next experiment.

### Pool depth is real and small

`eb8946e` multiplies the `slack` term by 16, the only term that adds blocks nobody has a claim on:
16808 -> 88808 blocks, 1.10 -> 5.83 GB, and 315 M points in flight at `--block 4096` against the
275 M that `--block 16384` has.  Interleaved against `63be9c4`, twice each:

| `--block 4096`, 150 + 105 | pool | routed | ns/point |
| --- | --- | --- | --- |
| baseline | 16808 blocks, 1.10 GB | 8.2 / 8.2 G | 18.3 / 18.3 |
| `slack * 16` | 88808 blocks, 5.83 GB | **9.0 / 9.0 G** | 16.6 / 16.7 |

**+10 %, reproducible, and it costs 4.7 GB.**  So depth in points is a real effect and a bad trade,
and it is not the mechanism: it closes a sixth of the 8.2 -> 13.0 G gap.  The branch stays unmerged.

### The hypothesis that fits the arithmetic

What survives is the receiver's side, and it fits quantitatively.  Every block a receiver takes is a
**fresh sequential stream**: a new 65600-byte region, written by a sender with non-temporal stores so
it is in DRAM and in nobody's cache, on pages the receiver has not walked.  A stream costs a roughly
fixed number of misses to establish -- prefetcher warm-up plus the page walk -- and then runs.  That
fixed cost is amortised over `block_points` points, so a 4x bigger block pays it 4x less often per
point.  The measured extra is **0.019 cache misses per point** at 4096 over 16384, i.e. **104 per
block transition** -- and 104 lines is 6.6 KB, the right order for a prefetcher to lock onto a new
stream and for the TLB to fault in a block spanning several pages.  It also explains why the cores
and the memory system are idle *together* at 262 GB/s: a latency-bound stream restart limits
memory-level parallelism without consuming bandwidth.

**This is a hypothesis and is labelled as one.**  The test that would confirm or kill it is a miss
census by role: `perf stat -e cache-misses,dTLB-load-misses -C <receiver cpus>` against
`-C <sender cpus>` (the placer pins, and prints the layout, so the two sets are known), which says
directly whether the extra misses are on the side that reads blocks.  If they are, the levers are
the block's page footprint (huge pages for the pool, already partly there since `b803bc1`) and a
receiver-side prefetch of the next block's head while it drains the current one -- which is the
`receiver prefetch` idea measured at ~5 % and dropped on the laptop, and which this says should be
retried here, on a machine where the effect exists.

### Where the extra misses are: the recycling path, and they are misses, not contention

The receiver-stream hypothesis above is **wrong**, and a miss census kills it:
`perf record -a -e cache-misses -c 5000`, `--sort srcline`, at both block sizes, 8 s of a steady
round.  The misses are overwhelmingly on the *sender* side, and the two biggest sources cost
**exactly the same per point** at both block sizes -- they are the router's baseline, not the effect:

| source | share at 4096 | per point | share at 16384 | per point |
| --- | --- | --- | --- | --- |
| `common.hpp:90` the NT store's source loads | 34.0 % | 0.0238 | 40.9 % | 0.0234 |
| `workers.hpp:161` the private write-combining line | 22.0 % | 0.0154 | 28.3 % | 0.0160 |
| `workers.hpp:164` | 12.0 % | 0.0084 | 13.8 % | 0.0079 |
| `workers.hpp:128` the destination word | 4.2 % | 0.0029 | 4.1 % | 0.0023 |

The whole difference is in the tail, and the tail is the **per-block recycling path**:

| source | 4096 | 16384 |
| --- | --- | --- |
| `workers.hpp:44/49/54/60` -- `free_pop_many`'s walk of the free ring's cells | 4.92 % | under 0.25 % |
| `atomic_base.h:536` -- the ring's read-modify-write | 3.49 % | 0.37 % |
| `workers.hpp:242` -- `Router_Grab` | 3.49 % | 1.56 % |
| `workers.hpp:109` -- `seal`'s write of `blk_link` | 1.61 % | 0.47 % |
| `ring.hpp:48` -- a receiver's inbox ring | 1.18 % | 0.27 % |

About 0.009 misses per point, ~36 per block transition, in the free ring, the seal's link and the
inbox -- every one of them a **cold shared line**, last touched by a core in another L3 domain or
another socket.  That is the mechanism, and it is a *miss* mechanism, not a *contention* mechanism.
Which is exactly why the two interventions aimed at contention failed: `762c660` removed the
service's polling and `c256db5` removed the seal's CAS convoy, but neither changed how many cold
lines a block costs to recycle.  `eb8946e`'s deeper pool helps by 10 % for the same reason it is
only 10 %: more blocks means a sender spins less in `free_pop_many`'s retry, but each pop still
walks the same cold cells.

**So the standing "per-group free list" item is the experiment to run, and my group-seal branch
sharded the wrong structure.**  Give each group its own free ring, so the block ids a sender pops
were pushed by a receiver in its own cache domain and the cells it walks are warm.  The seal's
`blk_link` and the inbox ring are the same argument.  That is the next session's work, and unlike
everything in sections 2 and 6 it is predicted by a measurement of where the misses are rather than
by a story about where the time must be going.

## What to change

1. **Find the per-block cost.**  §2 measures it and §6 rules out the service thread, which was the
   obvious suspect and is not the answer.  What is left to price, per block rather than per point:
   the single-word contention on `sealed_top` (one CAS per block from any of 150 senders) and on
   the free ring's `free_in`/`free_out`, and the sender's stall in `refill()` when the pool runs
   dry.  All three are node-wide singletons, which is what FINDINGS' standing "per-group free
   list" item was always about.  The next measurement is the sender's own per-point cost broken
   down, not another guess.  **Unmeasured: this is the next experiment on this node.**
2. **Until then, raise `Router_Opts::block_points` from 4096 to 16384.**  It is worth 60 % on one
   node and 67 % on two, and it is one number.  The cost is the pool: 1.10 -> 4.41 GB for a
   255-worker team, which the engines must subtract from what they hand the dictionary.  If that
   is too much, **8192 is the frugal setting**: 98 % of the single-node gain (12.7 of 13.0 G) at
   half the memory, but only 88 % of the two-node one (2.2 of 2.5 G) -- across the network the
   message size still matters on its own, so 16384 is the right default where the RAM is there.
3. **`service N % busy` counts only MPI work** and printed 0 % for a thread running at 2.19 GHz
   with 46 % of its cycles in one function.  It should count turns that did anything, or be
   dropped: as it stands it actively hid this bottleneck across sessions 4, 6 and 11.
4. The pool's *block count* is fixed by the team's shape, so `block_points` is the only lever on
   its footprint.  If 4.41 GB is unacceptable, make `n_blocks` shrink as `block_points` grows
   (the slack term `32 * S` and `R * (inbox + 1)` are both counted in blocks and are really
   "points in flight" budgets).  Unmeasured.

## Reproducing

```bash
oarsub -I -p grdix --project cryptanalyse --resource host=2,walltime=2
module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0
export HWLOC_ROOT=/gnu/store/hl3isj962ghsvvxig6ls84mw54nfc917-hwloc-2.13.0-lib
cmake -S . -B build -DCMAKE_BUILD_TYPE=release && make -C build -j64 router_bench

# one node
mpirun -np 1 --bind-to none build/examples/router_bench \
    --senders 150 --receivers 105 --rounds 2 --seconds 10 --block 16384

# two nodes, over HDR200
mpirun -np 2 --npernode 1 -H grdix-8,grdix-9 --bind-to none -x LD_LIBRARY_PATH -x PATH \
    --mca pml ucx -x UCX_NET_DEVICES=ibp66s0:1 --mca btl ^openib --mca oob_tcp_if_include br0 \
    build/examples/router_bench --senders 150 --receivers 105 --rounds 2 --seconds 10 --block 16384

# the fabric's own rate (server on grdix-9: ucx_perftest)
UCX_NET_DEVICES=ibp66s0:1 ucx_perftest grdix-9 -t tag_bw -s 262144 -n 20000 -w 1000

# the service thread (it is pinned to CPU 0), while a round is running
sudo-g5k sysctl -w kernel.perf_event_paranoid=-1 kernel.kptr_restrict=0
perf record -C 0 -F 999 -g -o svc.data -- sleep 12
perf report -i svc.data --stdio --no-children -g none --sort srcline
```

`--mca btl ^openib` only silences an "error initializing an OpenFabrics device" warning; the ucx
PML does not use that BTL.  Without `-x LD_LIBRARY_PATH -x PATH` the remote `orted` cannot find
the Guix Open MPI.

# Session 15 -- the direct engine at n = 40 on grdix: the dictionary is the wall, the push is 21 cycles

2026-09-13, **grdix-5** then **grdix-13** (2 x EPYC 9754, Zen 4c, 256 cores / 512 PU, 2 NUMA nodes,
32 L3 of 16 MB, 1003 GB), `oarsub -p grdix --project cryptanalyse -l nodes=1,walltime=2` (jobs
6924657 and 6924670).  Commit `86d69ca` on `bench-grdix` = `04d1041` of `router-fat-slack` (the
16x pool slack EXPERIMENT, `block_points` 16384, `is_good_pair` before `f(x)`) plus session 14's
FINDINGS; release build (`-O3 -DNDEBUG -march=native -ggdb`), debian13, `module load
openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0`.

**This session measures the application, not the router alone**: `double_speck64_demo --n 40
--ram 960G --producers-per-node 63 --dicts-per-node 192`, one rank, killed 30 s into the first
PROBE phase (a round is 34 s of FILL and about 590 s of PROBE; 19 rounds make the search).  The
question was where the Router's time goes when the senders run Speck and the receivers probe a
960 GB dictionary, which is the Router's design goal ("minimum time in Push and Pop with work
around them"), and `router_bench` cannot pose it.  Three variants of the engine, built from a
scratch worktree and never committed, separate the costs: **A** receivers pop and discard, **B**
producers compute, hash and pick a destination but never push, **C** the dict thread reads its
block in place and prefetches the slot of the point D ahead.

**Two nodes were lost to the method.**  A SIGKILLed 960G demo gives its memory back in about 5 s
(`MemAvailable` 849 -> 993 GB, the process a zombie meanwhile, logged twice on grdix-2), but once,
after the D = 16 run on grdix-13, **908 GB were still held two minutes after the kill** -- not
explained, and the reason a launch guard has to block rather than time out.  On grdix-5 a run
script's `pgrep -x double_speck64_` caught the previous run's dying process, declared its own run
dead without killing it, and the next run allocated a second 960G.  On grdix-13 the guarded script
waited 120 s for `MemAvailable`
after killing the D = 16 run, saw **85 GB**, and launched the D = 32 run anyway -- the wait fell
through instead of aborting.  Both nodes went down at once and OAR terminated the jobs (6924657,
6924670).  `run2.sh` now finds the demo as mpirun's child (`pgrep -P`), kills it with SIGKILL on
every exit path, logs `MemAvailable` after the kill and refuses to launch until it is back above
985 GB.  Every number below was reproduced on grdix-13 before it died, the baseline agreeing with
grdix-5 to the digit; sections 4-5's last measurements ran on **grdix-2** (job 6924675).

## Summary

| | FILL (6.0e10 inserts) | PROBE, steady | producer, cycles/point |
| --- | --- | --- | --- |
| **baseline** (`04d1041`) | 34.6 s, **1.8 G/s** | **1.8 G points/s** | 108 (throttled) |
| A: receivers discard | 21.4 s, 2.8 G/s | **3.5 G/s** | 56 = work + push |
| B: producers never push | 10.9 s, 5.5 G/s | **5.5 G/s** | 35 = work alone |
| C: prefetch D = 8 | 34.0 s, 1.8 G/s | **2.4 G/s** | -- |
| C: prefetch D = 16 | 33.6 s, 1.8 G/s | 2.4 G/s | -- |
| C: prefetch D = 32 | lost with grdix-13 | -- | -- |

The machine runs at the receivers' pace: 1.8 G points/s is 192 dict threads at one DRAM miss per
probe, 330 cycles a point at IPC 0.43, and the producers, able to push 3.5 G/s into idle
receivers, spend the difference throttled by the Router's backpressure.  The Router's own cost on
a producer is **21 cycles a point** (56 - 35), a third of the producer's time once the receivers
are out of the way -- not the 12 % the profile attributed to the push, and not nothing.  On a
receiver it is bounded above by the pure bench's pop cost, ~25 cycles of the 330.

## 1. The run

960 GB of dictionary is 120 G slots of 8 bytes in 192 shards of 5.0 GB, fill 0.5, so a round
inserts 6.0e10 = 2^35.8 preimages and probes all 2^40; 19 rounds.  The shards came up on
transparent 2 MB pages (`MAP_HUGETLB` refused: `nr_hugepages` = 0).  The Router built 47135
blocks of 16384 points (12.4 GB: the 16x slack is 32256 of them), private lines of 256 points
(0.8 MB per sender, the L2 rule), inboxes of 64 blocks, 32 groups of 1-2 senders and 6 receivers,
one per L3, 96 receivers and 31-32 senders per NUMA node, the service on CPU 0.

The live line settles at **1.8 G points/s in both phases**.  FILL takes 34.4-34.6 s, and its live
line shows the first 5 s doing almost nothing (136 M/s cumulative at 5.7 s, 730 M/s at 7.7 s) and
~2.0 G/s afterwards: the direct engine has no startup barrier, the dict threads are zero-filling
5 GB each when the round opens, and the producers fill 192 inboxes and drain the pool waiting for
them.  Round 0 only; 5 s of a 620 s round.

FILL reports **691277 blocks "held back for a full destination"** of the 3.66 M it moves: the
service could not place one block in five at its first try because the destination's inbox was
full.  That is the receiver-bound regime: the pool is in the inboxes and the parked lists, and the
senders wait in `refill()` for a receiver to release something.

## 2. Where the cycles go (perf, 20 s of PROBE, 199 Hz, every CPU)

The placer pins, so `-C` is an exact role filter (`/proc/<pid>/task/*/status` gives each thread's
CPU; the omp thread order is service, 192 receivers, 63 producers).  Counters over the 20 s:

| role | cycles | IPC | GHz | cycles per point |
| --- | --- | --- | --- | --- |
| 192 dict threads | 11.9e12 | 0.43 | 3.10 | **330** (142 instructions) |
| 63 producers | 3.90e12 | 0.97 | 3.10 | **108** |
| service | 6.19e10 | 1.10 | 3.10 | -- (110 K blocks/s) |

**Dict thread** (99 % in `dict_round`): `workers.hpp:236`, the first line of `Router_Pop`, 36.7 %;
`dict.hpp:156`, the tag compare on `A[s]`, 25.6 %; `:169-170` the counters 8.4 %; `:155` the
load 5.2 %; and **12.7 % in scalar Speck** (`double_speck64_problem.hpp:20-31, 143`): with m = 40
and 6.0e10 entries, **5.5 % of probes hit a genuine collision**, each one a key schedule and an
encryption in `is_good_pair`.  The pile on the head of `Router_Pop` is skid: it is the next
instruction to retire after the probe's stalled load, on L1-resident state that cannot cost 120
cycles by itself.  One DRAM miss per point, nothing in flight behind it: 330 cycles is the loaded
latency.  `cache-misses` reads 2.8 per point (whatever Zen 4's generic event counts).

**Producer** (99.5 % in `producer_round`): vector Speck (`double_speck64_problem.hpp:48-65, 137`,
the intrinsics) 44 %; `murmur64` of the image plus the destination pick and the lane loop
(`tools.hpp:94-96`, `producer.hpp:47-51`) 17 %; the push proper -- the line write
(`workers.hpp:158-164`) 7.4 %, `router_stream_copy` 2.5 %, the `fetch_add` on the destination
word 1.6 % -- **11.5 %**; **waiting for a free block** (`refill`, `workers.hpp:86`, and the
`pause` in `cpu_relax`) **19.7 %**.  Section 3 shows the 44 % of Speck is inflated: B does all of
it, the hashing included, in 35 cycles.

**Service** (90 % in `Router_Progress`, no MPI peer): the `retry_parked` loop over 193 targets
(`service.hpp:426-427, 129, 133`) 38 %; the inbox `push` that fails on a full inbox
(`ring.hpp:48-53`, which loads the receiver's `head` line) 20 %; `place()` 12 %; `PMPI_Testsome`
and `ompi_request_default_test_some` 8 % with nobody to talk to.  It hands over 110 K blocks/s
against session 14's 2.0 M/s ceiling: not a bottleneck, a core spinning on full inboxes and
touching every receiver's inbox head once per turn.

## 3. Two interventions: the receivers set the pace, the push is 21 cycles

**A -- receivers pop and discard** (`ctr[N_PROBE] += 1` in place of the probe; the FILL inserts
are discarded too).  PROBE runs at **3.5 G points/s**, FILL in 21.4 s (2.8 G/s; ~3.6 G/s once the
5 s ramp is past) with 70 K blocks held back instead of 691 K.  This is the producers-plus-Router
ceiling: 55.5 M points/s per producer, **56 cycles a point**.

**B -- producers never push** (`sink ^= h ^ dest` in place of `Router_Push`, the hash and the
destination pick kept).  FILL in **10.9 s**, PROBE at **5.5 G/s** (13.8 % of 2^40 at 27.5 s):
87 M points/s per producer, **35 cycles a point** for the vector Speck, the murmur and the pick.
(On grdix-5 the same binary took 17.7 s, launched 3 s after A's 960G was killed and while the
kernel was still giving it back; grdix-13 ran it twice at 10.9 s.  Never overlap two runs.)

So at this operating point **`Router_Push` costs a producer 56 - 35 = 21 cycles a point**, 6.8 ns
at 3.1 GHz, **37 % of the producer's time when the receivers keep up**.  The profile's 11.5 % was
wrong by 3x: the push's own instructions are few, and what they cost shows up as stall on
whatever comes next, which is Speck.  Session 14's pure bench measured 11.7 ns a point on a
saturated sender; with real work in the loop the marginal cost is 60 % of that, not the 10-20 %
the "it overlaps with compute" argument hoped for.  **Where the 21 cycles go is not settled** --
hypothesis, section 5.

And the baseline's producer at 108 cycles is neither: it is a 56-cycle producer throttled to the
receivers' 1.8 G/s, the other 52 spent in `refill` and in stores that wait.  Which is why the
question "what could be shaved in the Router" has two very different answers for the two sides.

## 4. The receiver: one miss in flight, and a prefetch buys 33 % in PROBE and nothing in FILL

**C -- `dict_round` reads the block in place** (`Router_Grab` / `Router_Release`), and before
point i it hashes point i + D and `__builtin_prefetch`es its first slot; inserts get the
write-intent form (`prefetchw`).  The murmur is computed twice per point, once for the prefetch.

| D | FILL | PROBE |
| --- | --- | --- |
| none (baseline) | 34.6 s | 1.8 G/s |
| 8, read prefetch in both phases | 33.8 s | 2.4 G/s |
| 8, write-intent in FILL | 34.0 s | 2.4 G/s |
| 16 | 33.6 s | 2.4 G/s |
| 32 | (the run that took grdix-13 down) | -- |

PROBE goes from 330 to **248 cycles a point** and stops there: D = 16 and 32 buy nothing over 8,
so what remains is not misses waiting one behind the other.  FILL does not move at all, with or
without write intent, although the same prefetch reaches the same shard.  Two readings, neither
measured yet: (a) the inserts' *writes* are the limit -- 2 G random dirty lines/s is 128 GB/s of
random write-back on top of 128 GB/s of fills, and random DRAM traffic is far from STREAM's
442-561 GB/s; (b) the receivers' DRAM traffic and the producers' non-temporal stores share the
memory controllers, and the whole node sits at a random-access ceiling around 2.0-2.4 G lines/s.
`randread` (below) puts a number on the read side of that ceiling.

**`randread2`** -- 192 threads pinned one per core, a private 3 GB array each on 2 MB pages, one
8-byte word read per uniformly random 64-byte line, which is what a probe does.  With *independent*
addresses the machine delivers **7.6-7.8 G random line fetches/s** whatever the number of streams
per thread (the core overlaps them on its own), and the demand-fill counter says they are real:
30.6 G `ls_dmnd_fills_from_sys.dram_io_near` for the run, 7.4 G of them the first touch, **1.00
DRAM fill per load** over the 23.2 G loads of the timed part.  That is 500 GB/s of line traffic,
90 % of STREAM's 561 GB/s read, on this pair of 12-channel DDR5 sockets with the TLB out of the
way.  With *dependent* addresses -- a stream's next address is the word it just loaded, so a
stream holds exactly one miss in flight:

| misses in flight per thread | G lines/s, 192 threads | per thread | latency seen |
| --- | --- | --- | --- |
| **1** | **1.41** | 7.3 M/s | **137 ns** |
| 2 | 2.55 | 13.3 M/s | 151 ns |
| 4 | 4.65 | 24.2 M/s | 165 ns |
| 8 | 7.14 | 37.2 M/s | 215 ns |
| 16 | 7.63 | 39.8 M/s | 402 ns |

One miss in flight per thread is 1.41 G/s, and **the baseline dictionary at 1.8 G/s (9.4 M probes/s
a thread, 106 ns a probe) is that regime**, a little better because the core overlaps the tail of
one probe with the head of the next.  Eight in flight reach the machine.  So the prefetched
receiver's plateau at 2.4 G/s is **not** the memory system: with D = 8 it could have 7 G/s, and what
stops it at 248 cycles a point is in the dict thread's own instruction stream (the collision check's
scalar Speck, the run loop's data-dependent branches, the two murmurs, Grab/Pop).  The perf profile
of the C variant was not taken and is the next measurement.  For scale, the same loop in L3
(2 MB a thread) runs 62-79 G/s independent and 12.8 G/s dependent at 15 ns, in L2 (256 KB) 79-116
and 22.7 G/s at 8 ns.

## 5. The sender: 21 cycles a point, hypothesis and test

A producer's private lines are F = 192 buffers of 256 points, 4 KB each, 768 KB of its 1 MB L2,
and each `Router_Push` touches two cache lines of them: the 64 B chunk being written (4 points a
line, a new one every 4 pushes to that destination) and the line's count word, which lives in
the *last* point of the 4 KB buffer, 4 KB away.  The hot set is 192 x 2 = 384 lines = 24 KB of a
32 KB L1D, before Speck's stack arrays and round keys.  **Hypothesis:** the 21 cycles are one or
two L2 hits per push (Zen 4: ~14 cycles each) from L1 thrashing.  If so, `--swc` will not help
(the count line is one per destination whatever the line size) but moving the counts into a
compact array (192 x 4 B = 12 lines, always hot) would, and so would a line size that keeps the
active chunk of every destination in L1.

**Measured** -- variant A (receivers discarding), `perf stat` on the 63 producer CPUs over 10 s of
PROBE:

| `--swc` | private lines per sender | PROBE | IPC | L1D misses per point | served by L2 |
| --- | --- | --- | --- | --- | --- |
| 64 | 192 KB | 2.9 G/s | 1.49 | 1.6 | 98 % |
| **256** (auto, the L2 rule) | 768 KB | **3.5 G/s** | 1.77 | **1.63** | 98 % |
| 1024 | 3 MB | 3.3 G/s | 1.69 | 1.7 | 71 % |

1.63 L1D misses a point, 98 % of them L2 hits (~14 cycles each on Zen 4, partly overlapped), is
the 21 cycles.  The hypothesis survives its first test: the miss count does not move with the line
size, as it predicts, since the count word is one line per destination whatever the size; 64 loses
on 4x more `stage_line` calls, 1024 loses on L2 misses (3 MB of lines in a 1 MB L2), and the auto
rule picked the best of the three.  **Not tested: the counts in a compact array** (192 x 4 bytes,
12 lines always hot), which the hypothesis says removes about half the misses, and a smaller line
that keeps every destination's active 64-byte chunk in L1 next to them.

## 6. `--benchmark` is not the producers' ceiling

`double_speck64_demo --n 40 --benchmark --producers-per-node 63` reports **20.3 M f/s per
producer, 1.3 G/s** -- 4.3x under B's 87 M/s per producer.  `benchmark()` walks one dependent
chain per lane (x = f(x)), so it measures the *latency* of a vector evaluation; the engine's
producer enumerates independent inputs and the core overlaps them.  As the denominator of "engine
rate over f/g rate" it makes the direct engine look 140 % efficient.  Either the benchmark should
enumerate independent inputs the way `producer_round` does, or the B variant of this session (the
engine's own loop minus the push) is the honest ceiling.

## 7. The service thread with one node

Nothing to shave on the critical path: 110 K blocks/s is 5 % of what it can do.  But it spins
through 193 parked lists and fails 192 inbox pushes per turn, each failure an acquire load of a
line the receiver writes, and polls MPI with no peer, 8 % of its cycles.  Harmless at np 1;
at np > 1 the same loop runs beside the real MPI work, unmeasured here.

## What to change

1. **The dictionary thread, not the Router, is where the machine's time goes**: 330 cycles a
   point, one miss, nothing behind it.  The prefetch (C) is 20 lines and is worth 33 % of PROBE
   on this machine, i.e. 3 hours of the 3.3-hour n = 40 search become about 2.3.  Where it stops
   (248 cycles, D-independent) is not DRAM: the memory system gives 192 threads 7 G random fetches/s
at 8 misses in flight each (section 4), so it is the dict thread's own code.  Profile the C variant.
Do it
   properly: carry the hash from the prefetch to the probe instead of computing murmur twice,
   and keep `Router_Pop` for PCS.  The memory note said "retry the prefetch on the cluster"
   (laptop: ~5 %); on 192 threads it is not 5 %.
2. **The push costs 21 cycles a point with room to push into**, 37 % of a producer.  Test the L1
   hypothesis of section 5 (`--swc`, then counts in a compact array); the Router's own knob that
   moves the lines' footprint is `swc_linesize`, and the L2 rule that sizes it may be the wrong
   rule -- the L1 is what the push touches per point.
3. **FILL's inserts do not gain from the prefetch**, PROBE's probes do.  Find out why before
   building on (1): if it is the write-back traffic, the FILL phase is at the memory wall and
   only fewer or smaller writes help (a 4-byte slot? a generation bit instead of the flush?).
4. **A 5 s ramp at round 0**: no startup barrier, the shards' zero-fill overlaps the first
   pushes.  One `#pragma omp barrier` after the dict threads' constructors, or nothing: 5 s once.
5. **`--benchmark` underestimates the direct engine's producers by 4x** (section 6).  Change
   the benchmark or stop quoting it as the ceiling.
6. **The service thread at np 1**: a harmless spinner, but `retry_parked` could skip receivers
   with nothing parked (a bitmap) and the MPI poll could be skipped with one node.  Unmeasured
   and off the critical path.
7. **Method:** one 960G demo at a time.  `pgrep -P` for the demo, SIGKILL it, then *wait for
   `MemAvailable`* -- the release is normally 5 s and was once over two minutes -- and abort rather
   than launch on a timeout.  A perf record over 512 CPUs takes 20 s to write out after `sleep`
   returns; time the kill from the PROBE mark, not from perf's exit.

## Reproducing

```bash
oarsub -p grdix --project cryptanalyse -l nodes=1,walltime=2 "sleep 7200"   # then ssh the node
module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0
cmake -S . -B build -DCMAKE_BUILD_TYPE=release \
    -DHWLOC_INCLUDE_DIR=/gnu/store/hl3isj962ghsvvxig6ls84mw54nfc917-hwloc-2.13.0-lib/include \
    -DHWLOC_LIBRARY=/gnu/store/hl3isj962ghsvvxig6ls84mw54nfc917-hwloc-2.13.0-lib/lib/libhwloc.so
make -C build -j 96 double_speck64_demo
sudo-g5k sysctl -w kernel.perf_event_paranoid=-1 kernel.kptr_restrict=0

# the run: FILL ~34 s, kill ~30 s into PROBE; build/session15/run2.sh does this with the guards
mpirun -np 1 --bind-to none --mca btl ^openib build/examples/double_speck64_demo \
    --n 40 --ram 960G --producers-per-node 63 --dicts-per-node 192

# roles: the placer pins, so thread -> CPU from /proc, omp order = service, 192 dicts, 63 producers
for t in /proc/$(pgrep -P $MPIRUN_PID -x double_speck64_)/task/*; do
    echo "$(basename $t) $(awk '/Cpus_allowed_list/{print $2}' $t/status)"; done | sort -n | grep -vE '[-,]'
perf record -a -F 199 -o probe.data -- sleep 20        # during PROBE
perf report -i probe.data --stdio --no-children -C <cpus of one role> --sort srcline
perf stat -C <cpus> -e cycles,instructions,cache-misses,cache-references -- sleep 20

# the variants: ~/mitm-exp on the Nancy home, -DEXP_A / -DEXP_B / -DEXP_C (EXP_D=<lookahead> at run time)
cmake -S ~/mitm-exp -B ~/mitm-exp/build-C -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS=-DEXP_C ...
# the DRAM random-access ceiling: independent and dependent (one miss in flight per stream) random lines
gcc -O2 -fopenmp -march=native build/session15/randread2.c -o randread2
OMP_NUM_THREADS=192 OMP_PLACES=cores OMP_PROC_BIND=spread ./randread2 3145728 2 1,2,4,8,16 b
perf stat -a -e ls_dmnd_fills_from_sys.dram_io_near,ls_dmnd_fills_from_sys.dram_io_far -- ./randread2 3145728 3 8 i
```

`build/session15/` in the grdix worktree keeps every log, tids file, perf.data and script of this
session (`run.sh`, `run2.sh`, `report.sh`, `batch2.sh`, `batch3.sh`, `check.sh`, `randread.c`,
`randread2.c`).

# Session 16 -- `--prefetch` on grdix: the ring's lookahead sweep

2026-09-13, **grdix-2** (2 x EPYC 9754, 1003 GB), `oarsub -p grdix --project cryptanalyse -l
nodes=1,walltime=1:30` (job 6924742).  Commit `e88aa67` on `bench-grdix` = `51f6d0b` of
`router-fat-slack` -- the dict thread's ring of `--prefetch` points between `Router_Pop` and the
shard (session 15's variant C, rebuilt on `Router_Pop` with a ring that spans blocks and is drained
at round end) -- plus session 15's FINDINGS.  Release build, debian13, the modules of session 14.

Same run as session 15: `double_speck64_demo --n 40 --ram 960G --producers-per-node 63
--dicts-per-node 192`, one rank, killed 30 s into the first PROBE (`build/session15/run2.sh`, the
guarded version), one run per lookahead in the interleaved order 0, 8, 4, 16, 2, 32, 12, 64, 24, 48,
so that a drift over the 15 minutes would show as a zigzag rather than a trend.  FILL is the round
report's time for 6.0e10 inserts (it includes round 0's 5 s ramp, section 15.1); PROBE is the live
line's cumulative rate about 30 s in.  `--prefetch 0` is the old loop, byte for byte.

## The sweep

| `--prefetch` | FILL, 6.0e10 inserts | blocks held back in FILL | PROBE, from the completion % |
| --- | --- | --- | --- |
| 0 (the old loop) | 34.0 s | 656 K | **1.77 G points/s** |
| 2 | 33.9 s | 530 K | 2.08 |
| 4 | 33.5 s | 227 K | 2.34 |
| 8 | 33.8 s | 115 K | 2.33 |
| 12 | 34.0 s | 104 K | 2.36 |
| **16** | 33.8 s | 93 K | **2.39** |
| 24 | 33.7 s | 95 K | 2.39 |
| 32 | 33.5 s | 109 K | 2.37 |
| 48 | 33.5 s | 105 K | 2.34 |
| 64 | 33.8 s | 117 K | 2.35 |

The PROBE column is the live line's completed fraction of 2^40 over its elapsed time, resolved to
0.1 % of the domain, i.e. about 1.5 % of the rate at 30 s; the live line's own "G points/s" prints
to one decimal and cannot separate 8 from 16.

**PROBE gains 34 %** (1.77 -> 2.39 G points/s) and the curve is a knee, not a hill: 2 buys half of
it, 4 almost all, and **4 to 32 is flat within the 2 % the measurement resolves**.  48 and 64 lose
about 2 %, at the edge of significance.  **16 is the default from here** (`51f6d0b` shipped 8; the
follow-up commit moves it), the middle of the flat range and so the choice least likely to fall off
it on a machine with a longer loaded latency or a smaller L1.  The ring costs nothing visible
against session 15's in-place variant: 2.39 here against 2.4 there, one copy and a second murmur per
point.

**FILL does not move**, 33.5-34.0 s for every D, while the blocks the service held back for a full
inbox fall 7x (656 K -> 93 K): the receivers do keep up better in FILL, and FILL is still 34 s.  So
FILL's ceiling is not the insert's miss latency.  The candidates are the memory system -- an insert
dirties a random line, and at 2.0 G inserts/s that is 128 GB/s of random write-back on top of 128
GB/s of fills, where session 15's randread2 measured reads alone -- and the producers, which in FILL
run f where PROBE runs g.  Session 15's variant A (receivers discarding) did FILL in 21 s, so the
producers are not it.  The random-write measurement is below.

## Where the prefetched receiver's cycles go

At `--prefetch 16` (`perf stat` on the 192 dict CPUs over 5 s of PROBE at 2.39 G points/s): **251
cycles a point, 215 instructions at IPC 0.86**, 28 branches of which **1.15 mispredicted per
point**, 2.2 L1D load misses per point.  Against D = 0's 330 cycles and 142 instructions at IPC
0.43, the prefetch took some 200 cycles of exposed miss out and put some 70 instructions in: a
second murmur for the prefetch address, the ring's copy and swap, the retire call.  What is left is
the instruction stream and two things it drags.  The mispredicts: the run loop's exit and the tag
compare are decided by the slot just loaded, a coin toss at load 0.5, so about one mispredict a
point at ~20 cycles on Zen 4.  The L1D misses: the prefetched slot lands in L1 but is gone 16 points
later, so the probe hits L2 (~14 cycles), and the block stream costs a line every 4 points.  The
source-line profile again puts 37 % at the head of `Router_Pop`, and it is skid again with a
different cause: a mispredicted branch's penalty lands on the resolved path's first instructions,
the counters after the loop and the loop top.  The 5.5 % collision rate's scalar Speck in
`is_good_pair` is 12 % of the instructions.

So in PROBE the receiver is no longer waiting on DRAM, it is executing: the node's random-read
ceiling is 7.6 G fetches/s (session 15) and the dictionary uses a third of it.  The lever left is
instructions per point -- carry the hash from the prefetch to the probe instead of computing murmur
twice, retire in place instead of copying and swapping, and test a whole 64-byte line of 8 slots at
once instead of walking the run slot by slot behind a data-dependent exit -- perhaps down to 150
cycles, ~4 G/s.  At that point the **producers**, 3.5 G/s with the push (session 15), are the wall,
and the push's 21 cycles the thing to shave.

## FILL is at the random-write wall

Three measurements settle why FILL never moves:

1. **randread2 in write mode** (one 8-byte word written per random 64-byte line, 192 threads, 3 GB
   each): **3.35 G lines/s** whatever the number of streams, against 7.5-7.8 G/s for reads in the
   same binary a minute later.  A random write is a fill and a write-back, two DRAM line operations,
   so 3.35 G writes/s is 6.7 G operations/s: the read ceiling in disguise.
2. **Half the dict threads, the same FILL.**  96 dict threads (and 159 producers) at `--prefetch 16`
   fill in **36.7 s**, 192 in 33.8 s; at `--prefetch 0` the 96 need 45.1 s.  With the prefetch,
   halving the threads costs 9 %: the limit is shared, not per thread.
3. **The blocks held back fall 7x** (656 K -> 93 K) as D grows while FILL stays at 34 s: the
   receivers keep up better, and something else does not.

At 2.07 G inserts/s (6.0e10 in the 29 s after round 0's ramp) FILL dirties 2.07 G random lines/s,
4.1 G DRAM operations, plus the block traffic (16 bytes written by the producer's non-temporal store
and read once by the receiver, 0.5 G lines/s each way): about 5 G operations/s against a 6.7-7.8 G
ceiling, in a mix of reads and writes that costs the bus turnarounds a pure read stream does not.
**FILL is at the memory wall for random writes, and the prefetch cannot lift it**; it only makes the
wall reachable with fewer dict threads.  An insert costs the memory system what 2.3 probes do.

## What to change

1. **`--prefetch` defaults to 16, shipped** (`7c90c48`; `51f6d0b` had 8).  Within 2 % of the best
   from 4 to 32; 16 matches the outstanding read requests an x86 core keeps in flight, which is the
   user's reason to expect it to travel to other x86 machines -- re-measure on grvingt (Skylake, 12
   fill buffers a core) before trusting it there.  The gain: PROBE 1.77 -> 2.39 G points/s, 34 %,
   the n = 40 round from about 620 s to about 490.
2. **The receiver's next 100 cycles are instructions, not misses** (section 2): carry the hash,
   retire in place, scan a line of 8 slots at once.  Then measure; 4 G/s would put the producers on
   the critical path.
3. **FILL is memory-bound on random writes**, 2.3 probes' worth of DRAM work per insert (section 3).
   Nothing in the dict thread fixes it; only fewer dirty lines would -- a layout where several
   inserts share a line, or a full-line write that skips the read-for-ownership.  Unmeasured.  And
   with the prefetch **96 dict threads fill nearly as fast as 192**: in FILL half the receivers are
   idle by construction, so the team's shape should be chosen for PROBE and accepted for FILL.
4. **PROBE at 2.39 G/s is 68 % of what the producers can push into idle receivers** (3.5 G/s,
   session 15). The two sides are 1.5x apart; after (2) the producer's 21-cycle push is the wall.
5. **Method**: the live line's one-decimal rate cannot separate lookaheads on a plateau; the
   completed fraction of the domain over the elapsed time can (`summ.sh`, the percentage column).  A
   profile that pins 37 % on the first line of a function on L1-resident data is skid, and the
   counters -- IPC, branch-misses, L1D misses -- say what it is skid from.

## Reproducing

```bash
# on the node, after the module loads and the release build of session 15
for D in 0 8 4 16 2 32 12 64 24 48; do
    EXTRA="--prefetch $D" TAG=pf$D STOP=probe bash -l build/session15/run2.sh
done
build/session15/summ.sh pf0 pf2 pf4 pf8 pf12 pf16 pf24 pf32 pf48 pf64     # the table
EXTRA="--prefetch D" TAG=pfDp STOP=probe PERF=1 STAT=1 bash -l build/session15/run2.sh   # the profile
build/session15/report.sh pfDp
```

# Session 18 -- the team's shape on grdix: 159 producers and 352 dict threads on the 512 hardware threads, 2.97 G points/s

2026-09-14, **grdix-2** (2 x EPYC 9754, Zen 4c family 25 model 160, 256 cores / 512 PU, 2 NUMA nodes, 32 L3
of 16 MB, 1003 GB), `oarsub -p grdix --project cryptanalyse -l nodes=1,walltime=3:00` (job 6926914).  Commit
`a1232e4` on `bench-grdix` = `fe7ae49` of `omp_reboot` (the dict thread's ring now prefetches in FILL as well as
in PROBE, `--prefetch` defaults to 8, the dict is a plain container whose `probe` fills a caller-owned vector)
plus sessions 14-16's FINDINGS.  Release build, debian13, `module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0`,
`kernel.numa_balancing = 0`, `perf_event_paranoid = -1`.

**The question**: with `--prefetch 8`, which split of the team between producers and dict threads runs the
direct engine fastest, on the 256 cores and on all 512 hardware threads, and where the best team's time goes.
Sessions 15 and 16 measured one split only, the 63/192 of the reference benchmark.

**The run**: `double_speck64_demo --n 40 --ram 900G --producers-per-node P --dicts-per-node R --prefetch 8`,
one rank, killed 30 s (the sweep) or 60 s (the finalists) into the first PROBE, one run per split, the 256- and
512-thread teams interleaved so that a drift would show as a zigzag (`build/session18/run18.sh`, `batch18.sh`).
**900G instead of the reference's 960G**: the Router's pool grows with the senders (its slack term is 16 x 32 x S
blocks of 256 KB, 41 GB at 255 senders) and 960G + 41 GB does not fit beside the OS in 1003 GB; per point
nothing changes -- the 63/192 anchor reproduces session 16 (PROBE 2.35 vs 2.33 G/s at `--prefetch 8`, FILL 2.06
vs 2.07 G inserts/s).  FILL is 5.6e10 inserts a round here (2^35.71), a round of PROBE 2^40 probes.

**Metric**: the live line's completed fraction of the domain over the elapsed time (0.1 % of 2^40 = 1.1e9 points,
so +-0.037 G/s at 30 s and +-0.018 at 60 s); the live line's own one-decimal rate cannot separate the splits, and
the 5-s and 10-s interval rates it allows move in steps of 0.22 and 0.11 G/s -- `summ18.py`'s `PROBElast` column
is that quantised number and misled for a moment.  FILL's rate is taken between the live line at 10 s and its
last, which leaves out round 0's zero-fill ramp (shorter with more shards).

## Summary

| team (P producers / R dict threads) | threads | FILL, G inserts/s | PROBE, G points/s |
| --- | --- | --- | --- |
| 63/192, the reference so far | 256 | 2.06 | 2.35 |
| 79/176 | 256 | 2.15 | 2.57 |
| **95/160**, the best on the cores alone | 256 | 2.16 | **2.61** |
| 111/144 | 256 | 2.10 | 2.56 |
| 127/384 | 512 | 2.15 | 2.87 |
| 143/368 | 512 | 2.13 | 2.91 |
| **159/352**, the best | 512 | 2.12 | **2.97** |
| 175/336 | 512 | 2.09 | 2.95 |
| 191/320 | 512 | 2.06 | 2.93 |
| 255/256 | 512 | 1.98 | 2.81 |

**The best team is 159 producers and 352 dict threads on all 512 hardware threads: 2.97 G points/s in PROBE, 26 %
over the 63/192 reference at the same prefetch, and 14 % over the best team on the cores alone (95/160, 2.61).**
The 512-thread plateau is broad: 143 to 191 producers within 2 %.  FILL is 2.1 G inserts/s for every team from
79/176 to 191/320, SMT or not: it is at the random-write wall of session 16, and both roles spin there.

The n = 40 round: FILL 31.5 s + PROBE 2^40 / 2.97e9 = 370 s = **6.7 min per round** at 159/352, against about
8.5 min at 63/192 with the same prefetch (the round report of `s159_352P` read PROBE 385 s, the stats' perf
runs included).

## 1. The sweep

The full sweep, 30 s of PROBE each (`s*` logs), then the finalists at 60 s (`s*L` logs):

| P / R | threads | FILL s | FILL G/s | blocks held back | PROBE 30 s | PROBE 60 s |
| --- | --- | --- | --- | --- | --- | --- |
| 31 / 224 | 256 | 43.3 | 1.46 | 38 K | 1.50 | |
| 47 / 208 | 256 | 35.0 | 1.86 | 66 K | 2.03 | |
| 63 / 192 | 256 | 31.8 | 2.06 | 117 K | 2.35 | |
| 79 / 176 | 256 | 30.0 | 2.15 | 486 K | 2.57 | |
| 95 / 160 | 256 | 30.0 | 2.15 | 2.0 M | 2.63 | 2.61 |
| 111 / 144 | 256 | 29.8 | 2.10 | 2.1 M | | 2.56 |
| 127 / 128 | 256 | 31.8 | 1.96 | 3.4 M | 2.32 | |
| 63 / 448 | 512 | 36.2 | 1.92 | 104 K | 2.23 | |
| 95 / 416 | 512 | 32.7 | 2.16 | 169 K | 2.76 | |
| 127 / 384 | 512 | 32.2 | 2.15 | 321 K | 2.87 | |
| 143 / 368 | 512 | 31.5 | 2.13 | 878 K | | 2.91 |
| 159 / 352 | 512 | 31.8 | 2.10 | 1.9 M | 2.99 | 2.97 |
| 175 / 336 | 512 | 31.3 | 2.09 | 1.9 M | | 2.95 |
| 191 / 320 | 512 | 31.2 | 2.06 | 2.1 M | 2.96 | 2.93 |
| 255 / 256 | 512 | 31.7 | 1.98 | 3.0 M | 2.81 | |

The "blocks held back" column is the Router's count of blocks the service could not place at its first try
because the destination's inbox was full -- the receiver-bound signature.  It grows a hundredfold from 31/224 to
127/128: the fewer the dict threads, the more the producers wait on them.  With `--prefetch`, 63/192 is on the
producer-bound side of the optimum (117 K held back, 37 M points/s per producer), not the receiver-bound side
session 15 found it on without the prefetch.

**On the cores alone the total is flat where the roles are balanced, and the dict thread's own rate is not.**
From 176 to 144 shards the dict thread goes from 14.6 to 16.3 to 17.8 M probes/s (sessions 15-16's per-thread
rate was taken as a constant), while the total stays at 2.56-2.61 G/s: the dict threads share their limit,
they do not each carry one.  The per-role counters below say what it is.

**On 512 threads the placer puts the receivers first**, one per core and then the second hardware thread of
each core, the senders after: the 159/352 team has 128 cores holding two dict threads, 127 holding a dict
thread and a producer, one holding a dict thread and the service (siblings are k and k + 256).  So SMT adds
dict threads to cores that already have one, and producers to the remaining cores' second thread; there is no
core with two producers.

## 2. Where the cycles go at the two best teams (perf stat, PROBE, 10 s per role)

Cycles and instructions are per hardware thread and per point of the machine's rate (`stat18.py`); the clock is
3.07-3.10 GHz on every role, SMT or not, so the cores are not throttled by the second thread.

| team | role | cycles / point | instructions / point | IPC | branch mispredicts / point | L1D misses / point |
| --- | --- | --- | --- | --- | --- | --- |
| 95/160 | dict thread (alone on its core) | **190** | 151 | 0.79 | 1.16 | 2.11 |
| 95/160 | producer (alone on its core) | 112 | 101 | 0.90 | | 1.95 (1.77 hit L2) |
| 159/352 | dict thread (SMT) | **367** | 150 | 0.41 | 1.16 | 1.90 |
| 159/352 | producer (SMT, next to a dict thread) | 164 | 99 | 0.60 | | 2.01 (1.63 hit L2) |
| 111/144 | dict thread | 174 | 150 | 0.86 | 1.16 | 2.11 |
| 191/320 | dict thread (SMT) | 339 | 153 | 0.45 | 1.16 | 1.90 |

Three things:

1. **The producer waits.**  Its work is 35 cycles a point and its push 21 (session 15), and it spends 112 to 201
   cycles a point depending on how many producers there are, at 100 instructions whatever the team: the rest is
   the push's wait for a free block.  Every team of the plateau is receiver-bound.
2. **A dict core delivers the same probes per second with one thread or two.**  190 cycles a point alone; two
   SMT dict threads at 367 cycles each are 184 core cycles a point.  The per-CPU instruction counts
   (`percpu18.py`, `perf stat -A`) show every dict thread of the 512 team at 8.3-8.5 M probes/s whether its
   sibling is a dict thread or a producer, and every dict thread of the 256 team at 16.4.  SMT buys the *machine*
   14 % because 352 dict threads at 8.4 M/s are more than 160 at 16.4, not because a core does more.
3. **The dict thread's instruction count fell from 215 (session 16) to 150** with the plain-container rewrite
   (`b6fb95d`, `fe7ae49`); at IPC 0.79 those are 190 cycles, of which the 1.16 mispredicts (the run loop's
   data-dependent exit at load 0.5) are some 25 and the exposed misses the rest.

## 3. What a probe costs the memory (the fill counters)

Zen 4's `ls_*_fills_from_sys` counters split the L1D fills by who asked -- a demand load, the hardware
prefetchers, a software `prefetch` -- and by where the line came from (`perf list`: `ls_dmnd_fills_from_sys`,
`ls_hw_pf_dc_fills`, `ls_sw_pf_dc_fills`, `ls_any_fills_from_sys`, each `.all / .local_l2 / .local_ccx /
.dram_io_near / .dram_io_far`).  On the 160 dict threads of 95/160, per probe:

| asked by | L1D fills / probe | of which from DRAM |
| --- | --- | --- |
| software prefetch (the run's first line, `dict.prefetch`) | 0.79 | 0.79 |
| demand loads | 0.25 | 0.23 |
| hardware prefetchers | 0.64 | 0.48 |
| **total** | **1.68** | **1.51** (1.49 near, 0.02 far) |

**A probe costs 1.5 DRAM lines**, not the 2.5-2.8 "cache-misses" of session 15 (a generic event that counts
something else).  The software prefetch and the demand DRAM fills add up to 1.0: the run's first line, brought
by the prefetch four times in five and by the load itself once -- the once is the run that crosses into a
second line (at load 0.5 a run is 2-3 slots long, so about a fifth of the probes read a second line the
prefetch did not ask for) plus the ring's lookahead falling short.  The block stream is 16 bytes a point, a
quarter line, and the hardware prefetcher fetches it -- it is sequential -- and about another quarter line of
neighbours nobody reads.  The far-DRAM share is 1 %: the placer's NUMA interleave and the shard's first touch
work.

Machine-wide (`perf stat -a`): **5.3 G DRAM lines/s in PROBE at 95/160** (2.02 lines per point over every CPU,
the producers' and the service's share included), against the 7.6-7.8 G/s random-read ceiling of session 15's
randread2; at 159/352 about 6.0 G/s of fills plus 0.75 G/s of non-temporal block writes.  **The 512-thread
PROBE runs the DRAM at 80-90 % of its random-access ceiling**; the 256-thread one at 70 %.  That is the shared
limit of section 1: the loaded latency at that utilisation is what the dict thread's eight outstanding
prefetches divide, and more dict threads (or a second per core) buy throughput only as far as the DRAM has
left.

## 4. FILL is at the random-write wall, and both roles spin there

Every team from 79/176 to 191/320 fills at 2.06-2.16 G inserts/s.  The counters at 159/352 during FILL: the
dict thread runs **210 instructions a point where FILL's insert path is 99** (95/160's FILL count), the producer
141 where its work is 91 -- both are spinning, the dict thread in `Router_Pop` for a block, the producer in
`Router_Push` for a free one, and neither is the limit.  Per insert the dict thread fills 2.66 L1D lines
(0.64 demand, 0.70 hardware prefetch, 1.32 software prefetch) of which **2.5 from DRAM**, and the machine
moves 5.9 G DRAM lines/s of fills plus the 2.1 G/s of dirty lines it writes back and 0.5 G/s of block writes:
about 8.5 G line operations/s, the 6.7-7.8 G ceiling of session 16 with the mix's turnarounds.  Nothing in the
team's shape moves it; the layout of the dictionary (several inserts per dirty line, or a full-line write) is
still the only lever, unmeasured.

## 5. The best team's profile (perf record, 20 s of PROBE, 199 Hz, every CPU)

`perf record -a -F 199 -- sleep 20` on `s159_352P`, 2.0 M samples, the roles separated by the CPU lists of
the pinned threads (`-C`), `--sort srcline` (`report18.sh`, on the node: the report over a 107 MB `perf.data`
takes minutes and does not belong on the frontend).

**Dict thread (352 SMT threads, 8.4 M probes/s each).**  94.3 % in the engine's `run` loop (everything inlined),
5.0 % in `Speck64128KeySchedule`.  By source line:

| share | line | what |
| --- | --- | --- |
| 25.7 % | `dict.hpp:188`, `ctr[N_PROBE] += 1` after `probe` | skid: the first instructions after the run loop's mispredicted exit |
| 19.1 % | `dict.hpp:145`, the tag compare on `A[i]` | the slot load's exposed part: the run's second line, the prefetch that fell short |
| 11.5 % | `dict.hpp:189`, `for (u64 x : candidates)` | skid again, the loop over an almost always empty vector |
| 3.9 % + 2.5 % + 1.2 % + 1.9 % | `dict.hpp:144, 142, 143, 148` | the `while (A[i] != 0)`, the tag, the index, the wrap |
| 13 % | `Speck64128KeySchedule` + `double_speck64_problem.hpp:20-31` | `is_good_pair`'s scalar Speck, on the 5.5 % of probes that have a candidate |
| 3.3 % + 0.9 % + 0.7 % | `workers.hpp:236, 242, 241`, `Router_Pop` | the pop, no longer a pile (37 % in session 16: that was the unprefetched miss's skid) |
| 2.7 % + 1.6 % + 1.1 % + 1.0 % | `dict.hpp:229, 234, 232, 224` | the ring |
| 1.3 % | `dict.hpp:215`, `cpu_relax()` | the spin on an empty inbox: a little starvation even here |
| 0.9 % | `tools.hpp:92` | the second murmur |

Two thirds of the dict thread is the probe loop and what its data-dependent exit costs; the levers of session 16
stand (test a line of 8 slots at once, no per-slot branch; carry the hash; retire in place), and the run's
second line is the one prefetch the ring does not issue.  The Speck at 13 % is the price of m = 40 with a
900 GB dictionary: 5.5 % of the probes hit a genuine 40-bit collision and each one is a key schedule and an
encryption -- `vfg` on a batch of candidates, or more check bits, would take most of it back.

**The same dict thread alone on its core (`s95_160P`, 16.4 M probes/s)** has the same shape with the
proportions shifted towards the compute: the probe loop and its skid (`:188, 145, 189, 144, 142, 143, 148`)
56 % instead of 66 %, `is_good_pair`'s Speck **28 %** instead of 13 % (`Speck64128KeySchedule` 10.7 % alone),
`Router_Pop` 3.3 %, the ring 7 %.  In cycles the Speck is the same 50 a point both ways (28 % of 190, 13 % of
367): the SMT thread's extra cycles are all stall, on the probe's lines, which is what a busier DRAM looks like
from inside the core.

**Producer (159 SMT threads, next to a dict thread each, 18.7 M points/s).**  99.5 % in `run`.  The vector
Speck (`double_speck64_problem.hpp:48-65, 137` and the AVX-512 intrinsics) about 35 %; the push -- the private
line's write at `workers.hpp:161` (13.0 %, the L1D miss served by L2 of session 15), its neighbours
`:158-164` (13 %) and the non-temporal stream copy of a full line into the block at `common.hpp:93` (6.0 %) --
about 32 %; `murmur64` and the destination pick (`tools.hpp:92-96`, `producer.hpp:39-43`) about 20 %.  The
producer's 164 cycles a point against its 56 of work do not show up as a spin: the waiting is the SMT sibling's
share of the core (IPC 0.60 here, 0.90 alone at 95/160) and the memory stalls of the line write and the
stream copy under a DRAM at 80-90 % load.  Alone on its core (`s95_160P`) the producer is half vector Speck
(50 %), a fifth push (22 %: the line write 5.8 %, the stream copy 3.8 %, and **6.5 % in the free ring -- the CAS
of `free_pop_many` at `workers.hpp:65` and `refill`'s spin at `:86`**, the wait for a block the receivers have
not released), a fifth murmur and lane loop (22 %).

**Service thread**: 81 % in `Router_Progress` (its own loop: `service.hpp:426`, the ring reads, the atomics),
16 % in `MPI_Testsome` -- a spinner with one node, 1.0-1.2 cycles per thousand points routed, 5 % of a core.

## 6. Switching the hardware prefetchers off (PrefetchControl, MSR C000_0108)

Zen 4 exposes a per-core **PrefetchControl** MSR (`C000_0108`; bit 0 L1 stream, 1 L1 stride, 2 L1 region, 3 L2
stream, 5 up/down; 1 = off; it reads 0 on the debian13 image, everything on).  `pf18.sh MASK` writes it on the
512 CPUs through `/dev/cpu/*/msr` under `sudo-g5k` and reads it back; `pfbatch18.sh` runs the best team under a
mask and restores 0 at exit.  Each run is 60 s of counters inside PROBE (about 110 s of it in all, +-0.01 G/s),
FILL counters included:

| PrefetchControl | PROBE, G points/s | vs default | dict thread, DRAM lines / probe | hw-prefetch L1D fills / probe | machine DRAM, G lines/s | FILL, G inserts/s |
| --- | --- | --- | --- | --- | --- | --- |
| 0x0, the defaults | 2.95 | | 1.46 | 0.40 | 5.90 | 2.10 |
| **0x2f, everything off** | **3.06** | **+3.7 %** | 1.35 | 0 | 5.77 | 2.13 |
| 0x4, L1 region off | 2.99 | +1.4 % | 1.38 | 0.21 | 5.71 | 2.11 |
| 0x8, L2 stream off | 2.97 | +0.7 % | 1.48 | 0.41 | 6.08 | 2.12 |
| 0x1, L1 stream off | 2.87 | **-2.7 %** | 1.48 | 0.42 | 5.85 | 2.10 |
| 0x20, up/down off | 2.95 | 0 | 1.46 | 0.40 | 5.92 | 2.11 |

**Everything off is 3.7 % faster.**  The L1 stream prefetcher earns its keep -- it is what reads the block
stream ahead of `Router_Pop`, and losing it alone costs 2.7 % -- the region prefetcher wastes a fifth of a line a
probe on neighbours of a random slot that nobody reads (+1.4 % without it), the L2 stream prefetcher's traffic
never reaches the L1D counters but shows in the machine's DRAM count, and with all of them off the block stream
comes in by demand (0.58 demand fills a probe instead of 0.37) and the machine still moves 2 % fewer lines for
3.7 % more points.  **That is the DRAM-bound signature: fewer lines per probe, proportionally more probes**,
which is the second confirmation, after section 3's 80-90 % of the random-access ceiling, that on 512 threads
PROBE is bound by the memory system and not by the core.  On FILL the prefetchers make no difference (2.10-2.13).
The MSR needs root and is per boot; it is a note for a BIOS setting, not a code lever.  The MSR was restored to 0
on every CPU at the end (`phase3.out`), and the FILL-phase "all CPUs" lines of `stat18.py` can straddle the
phase's end (`pf0x4`'s 3.38 GHz reading) and are not to be quoted.

## What to change

1. **The reference benchmark on grdix is now `--producers-per-node 159 --dicts-per-node 352 --prefetch 8` on
   all 512 hardware threads, `--ram 900G`** (the pool is 28.8 GB at 159 senders; 960G no longer fits), 2.97 G
   points/s in PROBE, 2.1 G inserts/s in FILL, 6.7 min per n = 40 round.  Without SMT, 95/160 (2.61).  The
   63/192 of sessions 15-16 is 26 % below the best and on the producer-bound side once the receivers prefetch.
   The plateau is broad (143-191 producers within 2 %), so the split need not be tuned to the unit.
2. **On 512 threads PROBE is DRAM-bound at 1.46 lines a probe**, 5.9 G lines/s of fills plus 0.75 G/s of block
   writes against a 7.6-7.8 G/s random-access ceiling, and the rate follows the lines per probe (section 6).  The
   lever is **fewer DRAM lines per probe**: a line-aligned layout of 8-slot buckets makes a probe exactly one line
   (no second line for the run that crosses: 0.2-0.3 demand lines a probe today), a compare of 8 slots at once
   takes the run loop's 1.16 mispredicts and part of the 150 instructions with it, and the same layout spares
   FILL its second line too.  Expected: 1.46 -> 1.25 lines a probe, up to +15 % PROBE if the DRAM stays the
   wall.  Session 17 asked for the same thing from grvingt's side ("8-slot buckets or the 8-slot compare would
   give fill 0.5 the fill 0.25 cost").
3. **On 256 threads PROBE is latency-bound** (70 % of the ceiling, the dict thread's 8 outstanding prefetches
   dividing a loaded latency), and there the dict thread's cycles matter: `is_good_pair`'s scalar Speck is 50
   cycles a point (28 % of the dict thread alone on its core) for 5.5 % of probes -- more check bits, or a
   vectorised batch of candidates, take most of it; the run's second line is the one prefetch the ring does not
   issue.  Both are moot on 512 threads until the lines per probe come down.
4. **The producers have slack at every team of the plateau**: 100 instructions a point, 56 cycles of work, 112 to
   200 cycles spent.  Session 15's 21-cycle push is not on grdix's critical path any more; leave it.
5. **FILL stays at 2.1 G inserts/s** for any team, SMT or prefetchers or not: 2.5 DRAM lines an insert plus the
   write-back, both roles spinning.  Only the dictionary's layout can move it (item 2 takes the second line;
   the dirty line per insert stays).
6. **The hardware prefetchers cost 3.7 %** on this access pattern (all off is fastest; keep L1 stream if only one
   is kept) -- a BIOS or per-boot root setting, not something the code does; noted, not recommended for the
   library.
7. **Method.**  (a) The live line's completed fraction over the elapsed time is the rate metric, at 60 s or more;
   its 5-s and 10-s differences are quantised to 0.22 and 0.11 G/s and mislead.  (b) An n = 40 PROBE lasts 370 s
   at this rate, and a `perf record -a` of 512 CPUs writing 107 MB to the NFS home took six minutes once
   (`s159_352P`), so the counters that followed it straddled the phase's end: write `perf.data` to local disk,
   or take the counters before the record.  (c) `perf report` over such a file runs on the node, never on the
   frontend.  (d) Zen 4's `ls_{dmnd,hw_pf,sw_pf,any}_fills_from_sys.*` are the counters that say what a probe
   costs the memory and who asked for it; the generic `cache-misses` (2.8 a point in session 15) counts something
   else.  (e) `perf stat -A` per CPU, joined with the pinned threads' CPU lists, separates the two siblings'
   roles.

## Reproducing

```bash
# on the node, after the module loads (build/session18/setup18.sh builds and sets the sysctls)
bash -l build/session18/batch18.sh                                    # the sweep, 12 splits, 30 s of PROBE each
SPLITS="159/352 95/160 191/320 111/144 175/336 143/368" TAGSUF=L PROBE_S=60 STAT=1 bash -l build/session18/batch18b.sh
SPLITS="159/352 95/160" TAGSUF=P PERF=1 STAT=1 FILLSTAT=1 PROBE_S=95 bash -l build/session18/batch18c.sh   # the profiles
MASKS="0x0 0x2f 0x4 0x8 0x1 0x20" bash -l build/session18/pfbatch18.sh   # the prefetcher intervention, restores 0
python3 build/session18/summ18.py s159_352L ...     # the rate table (light: fine on the frontend)
python3 build/session18/stat18.py s159_352P         # the per-role counters; percpu18.py TAG [fill] per thread
bash build/session18/report18.sh s159_352P          # perf report by symbol and source line per role -- ON THE NODE
```

# Session 20 -- four grdix nodes over HDR200: PROBE is the wire at 16 bytes a point, FILL loses a quarter to the pool

2026-09-15, **grdix-4, -5, -6, -7** (each 2 x EPYC 9754, 256 cores / 512 PU, 2 NUMA nodes, 1003 GB; one
Mellanox HDR200 port `ibp66s0`, 4X HDR = 200 Gb/s = 25 GB/s a direction, MTU 4096, fw 20.43.2566), job 6927307,
`oarsub -l {cluster='grdix'}/nodes=4,walltime=3:0:0 -q production`.  Commit `b13c3db` on `bench-grdix` =
`fe7ae49` of `omp_reboot` plus session 18's FINDINGS; the same release binary as session 18 (2026-09-14 22:40).
debian13, `module load openmpi/4.1.6 hwloc/2.13.0 ucx/1.20.0` (Guix Open MPI 4.1.6, UCX 1.20.0, transports
`rc_mlx5`/`dc_mlx5`), `perf_event_paranoid = -1` on grdix-5, the profiled node.

**The question**: the direct engine on several nodes for the first time on grdix.  Starting from session 18's
best team on the cores alone, 95 producers and 160 dict threads a node, what does a 4-node run deliver per node
against a node alone, and what binds it -- the wire, the service thread, the dict threads, the Router's flow
control?

**The run**: `double_speck64_demo --n 42 --ram 900G --producers-per-node 95 --dicts-per-node 160 --prefetch 8
--nrounds 1`, one rank per host, `--nrounds 1` so that the run ends by itself after round 0's PROBE with the
exact round report.  Session 18's problem-size hygiene keeps `fill * w / 2^m` at 5 %: four nodes hold four
times the slots, so n goes from 40 to 42 (FILL 2.25e11 inserts = 2^37.71, PROBE 2^42, a 10-minute round).  The
sweeps run at 225G a node with n = 38 / 39 / 40 for 1 / 2 / 4 nodes (5.1 % everywhere, a 3-minute round).
Open MPI over UCX on the HDR port:

```
mpirun -np 4 --npernode 1 -H grdix-4,grdix-5,grdix-6,grdix-7 --bind-to none -x LD_LIBRARY_PATH -x PATH \
    --mca pml ucx -x UCX_NET_DEVICES=ibp66s0:1 --mca btl ^openib --mca oob_tcp_if_include br0 ...
```

with `LD_LIBRARY_PATH` set to the Guix Open MPI's `lib` by hand: **the module leaves it empty and the demo has no
RUNPATH for `libmpi`, so without it the binary built against 4.1.6 loads debian13's `libopenmpi40` 5.0.7** (it
happens to run at np 1, which is how session 18 never saw it).  Scripts and logs in `build/session20/`: `run20.sh`
(one run over the first NP hosts of `nodes`, a memory-and-liveness guard on every host before it, a per-second
sample of every host's `port_xmit_data` / `port_rcv_data` / `port_xmit_wait` through `ibsample.sh`), `batch20.sh`
(the series), `summ20.py` (the exact round reports), `live20.py` (the plateau rates), `ib20.py` (the wire),
`roles20.sh` / `waitprof20.sh` / `waitprofF20.sh` (perf per role on one node, in PROBE or in FILL).

**Metric**: three, and they differ on purpose.  The **exact round report** gives the phase's time and its
2^k evaluations, ramp and drain included.  The **plateau** is the live line's completed fraction over 5-second
windows, the median over the middle 60 % of the phase (`live20.py`; the fraction has one decimal, so a 5 s window
of an n = 42 phase is quantised to +-11 % and the median over 60 windows is what one reads, +-3 % at n = 40).
The **wire** is the HDR port's own counters in 10-second bins, the plateau being the median of the bins above half
the peak: what actually crossed the link, IB headers included (about 1 % over the MPI payload, which the engine
prints as `node-->`).

## Summary

| run (225G a node unless said) | np | FILL exact, G inserts/s (a node) | FILL plateau a node | PROBE exact, G points/s (a node) | PROBE plateau a node | wire, GB/s a node, each way |
| --- | --- | --- | --- | --- | --- | --- |
| **900G, 95/160, n 42** | 4 | 5.10 (1.27) | 1.82 | **7.65 (1.91)** | **1.98** | **24.3** |
| 900G, 95/160, n 42, the first run | 4 | 5.4 (1.35) | -- | 8.0 (live line) | -- | 24.3 |
| 900G, 95/160, n 40, a node alone | 1 | 1.85 | 2.37 | 2.47 | 2.64 | -- |
| 95/160, n 38 | 1 | 1.85 | 2.32 | 2.52 | 2.58 | -- |
| 95/160, n 39 | 2 | 2.05 (1.03) | 1.07 | 4.16 (2.08) | 2.25 | 18.3 |
| 95/160, n 40 | 4 | 4.32 (1.08) | 1.53 | 7.53 (1.88) | 2.03 | 24.2 |
| 95/160, n 40, repeat | 4 | 4.04 (1.01) | 1.52 | 7.67 (1.92) | 2.03 | 24.5 |
| 63/192 | 4 | 4.43 (1.11) | 1.69 | 7.56 (1.89) | 2.03 | 24.5 |
| 159/352 (SMT) | 4 | 2.48 (0.62) | 0.60 | 6.94 (1.74) | 1.92 | 23.0 |
| 95/160 `--block 65536` (1 MB messages) | 4 | 3.58 (0.89) | 0.98 | 7.38 (1.85) | 2.03 | 23.3 |
| 95/160 `--block 4096` (64 KB messages) | 4 | 2.76 (0.69) | 1.01 | 3.69 (0.92) | 0.99 | 12.7 |
| 95/160 `--credit 16` | 4 | 3.51 (0.88) | 1.29 | 7.21 (1.80) | 1.92 | 23.8 |
| 95/160 `--n-recv 128` | 4 | 4.23 (1.06) | 1.53 | 7.51 (1.88) | 1.98 | 23.8 |
| 900G, 95/160, n 42, `--inbox 256` (pool 93663 blocks) | 4 | 5.02 (1.25) | 1.50 | 7.58 (1.89) | 1.98 | 24.3 |
| 900G, 95/160, n 42, `--n-recv 128 --credit 8` | 4 | 5.38 (1.34) | 1.81 | 7.68 (1.92) | 1.98 | 24.4 |
| 95/160, n 39, `--credit 16` | 2 | 1.97 (0.98) | 1.24 | 4.25 (2.12) | 2.31 | 18.4 |
| 95/160, n 39, `--block 65536` | 2 | 1.43 (0.72) | 0.69 | 4.20 (2.10) | 2.31 | 17.3 |
| 900G, 95/160, n 42, `--credit 2` | 4 | 5.40 (1.35) | 1.79 | 6.78 (1.69) | 1.76 | 21.4 |
| 900G, 95/160, n 42, **diagnostic binary: the senders leave 8192 blocks in the free ring** instead of 32 | 4 | 5.39 (1.35) | 1.81 | -- (killed after FILL) | -- | 24.3 |
| 900G, 95/160, n 42, `--inbox 8` | 4 | 5.42 (1.35) | 1.80 | 7.64 (1.91) | 1.98 | 24.1 |

**Four nodes run PROBE at 7.7 G points/s, 1.9-2.0 a node against 2.6 for a node alone (75 %), and the wire is
the reason: every node sends 24.3 GB/s and receives 24.3 GB/s for the whole phase, 97 % of the port's 25 GB/s.**
A point is 16 bytes on the wire and three quarters of them leave the node, so the port caps a node at
24.3e9 / (16 x 3/4) = 2.03 G points/s -- the plateau measured.  Nothing else moves it: 63/192 gives the same
24.5 GB/s, messages of 1 MB the same, more credit or more posted receives the same, and the 512-thread team is
slower, not faster.  **FILL runs at 1.8 G inserts/s a node against 2.4 alone (77 %), with the wire at 15 GB/s: it
is not the wire but the Router's pool**, which sits in sealed blocks waiting for the peers' credit while the dict
threads wait for points half of the time.

## 1. PROBE is the wire, at 16 bytes a point

The port counters of the four hosts over the 900G run's PROBE (580 s of plateau, 10-s bins): 24.3 GB/s out and
24.3 GB/s in on every host, bins within 0.2 GB/s, 13.6 TB moved each way per host, `port_xmit_wait` flat on
grdix-4/-5 and rising on grdix-6/-7 in the first run only.  HDR 4X is 200 Gb/s = 25.0 GB/s a direction: the
link is at 97 %, and the 1 % between the counters and the engine's 24.1 GB/s of MPI payload is the IB headers
at MTU 4096 and the rendezvous control traffic.  A link that full is the ceiling, and the arithmetic closes:

| | a node alone | 4 nodes | N nodes |
| --- | --- | --- | --- |
| points leaving the node | 0 | 3/4 | (N - 1) / N |
| wire cap, G points/s a node, at 16 B a point and 24.3 GB/s | -- | **2.03** | 24.3 / (16 (N - 1) / N): 1.62 at 16 nodes, **1.52 at N -> inf** |
| measured PROBE plateau a node | 2.64 | **1.98-2.03** | |

So a 64-node grdix job would run PROBE at 1.55 G points/s a node, 59 % of what the node does alone, whatever
the team, and **the only lever is the number of bytes a point costs on the wire**: 12 bytes lifts the cap to 2.7 a
node at 4 nodes (then the dict threads bind again at 2.6), 8 bytes to 4.05; at N -> inf, 2.0 and 3.0 against
today's 1.52.  The wire point is `Point {key, val}`, two u64: the key's top bits pick the destination and are
known to the receiver, the check bits are what the shard stores, and `val` is a preimage of n <= 42 bits.

**What does not move it** (all at 225G, 4 nodes, PROBE plateau a node / wire):
- **the team**: 95/160 and 63/192 both 2.03 / 24.5 GB/s -- the 4-node rate is team-independent where the
  1-node rate is not (session 18: 2.61 vs 2.35);
- **the message size**: `--block 65536` (1 MB messages, 66 GB of pool) 2.03 / 23.3;
- **the window**: `--credit 16` 1.92 / 23.8, `--n-recv 128` 1.98 / 23.8, within the run-to-run scatter (the
  two identical 95/160 runs differ by 2 %);
- **the 512-thread team**, 159/352: **1.92 / 23.0, 5 % under the 256-thread teams**.  The placer gives the
  service thread's core a dict thread as its SMT sibling (session 18), and the service thread is the one thread
  that must not lose cycles here: see 2.

**What halves it**: `--block 4096`, the 65600-byte message, 0.99 / 12.7 GB/s.  Four times the messages for the
same bytes, and the service thread stops at about 190 K messages/s each way: at 262 KB the same count would be
50 GB/s, so at the default block it has a factor two of slack and the wire binds first.

## 2. The service thread: zero-copy, 92 K messages a second each way, half loaded

`perf record -t` on grdix-5's main thread during PROBE (900G and 225G runs alike; IPC 1.3-1.5, 3.1 GHz):

| share | where |
| --- | --- |
| 35 % | `Router_Progress` itself: 17 % scanning the 163 parked lists (`retry_parked` over R + n_nodes targets, every turn), 5 % `complete()` over a block's 256 line counters, the rest the turn's own loop |
| 15 % | `MPI_Testsome` (`PMPI_Testsome` + `ompi_request_default_test_some`), two calls a turn over 32 slots each |
| 7.5 % | `uct_ib_mlx5_devx_mkey_pack` + 2 % `ucp_rkey_pack_memh` + 1.4 % `ucs_pgtable_lookup`: the memory key of each rendezvous, per message |
| 6.5 % | `ucp_tag_recv_nbx`, 5 % `ucp_tag_rndv_rts_progress`, 1.3 % `uct_rc_mlx5_base_ep_get_zcopy`: the receiver side of the rendezvous, an RDMA read per message |
| 0.13 % | `memmove` |

Every 262 KB message is one rendezvous: the sender posts an RTS, the receiver reads the block by RDMA into the
posted receive's block, and the service thread touches no byte of it.  Its cost is per message, and at 92 K
messages/s each way (24.3 GB/s / 262208 B) it is about half loaded, by the `--block 4096` run's 190 K/s wall.
Everything below 16 % of it is Open MPI's and UCX's own bookkeeping; the Router's own 35 % is mostly scanning
empty parked lists, a cost that grows with R and would be worth a non-empty list instead when the service thread
becomes the bound -- which it does at 65 KB messages, or on a faster wire, or with a dict thread on its sibling.

## 3. The dict threads and the producers wait for the wire

`perf record -C` over the 160 dict CPUs and the 95 producer CPUs of grdix-5, 8 s each, by source line (the
engine's functions are all inlined into `run`, so symbols say nothing):

| role, PROBE at 4 nodes | waiting | working |
| --- | --- | --- |
| dict thread | **55 %** (37 % `Router_Pop` finding the inbox empty, 18 % the `cpu_relax` behind it) | 45 %: the probe loop 19 %, `is_good_pair`'s Speck 9 %, the rest hashing and the ring |
| producer | **22 %** (15 % `refill`'s `cpu_relax` for a free block, 7 % the free ring's head) | 78 %: Speck 45 %, the push 10 %, the rest |

Both roles idle on the transport: the dict threads for want of points, the producers for want of blocks -- the
pool's blocks are sealed and parked for the peers' credit (194 M parkings in the 900G PROBE, 0 receives left
unposted), which is the healthy shape of a wire-bound run: the outbound queue is never empty and the inbound
side never lacks a block.

## 4. FILL: not the wire, and no transport knob moves it

FILL at 4 nodes runs at 1.82 G inserts/s a node on the plateau (1.27 exact, the 44 s phase carrying a ramp and a
drain of several seconds) against 2.37 alone: **77 %, with the wire at 15.3 GB/s, 61 % of the port**.  On
grdix-5 during the 900G FILL (roles at 30-46 s of the phase):

| role, FILL at 4 nodes | waiting | working |
| --- | --- | --- |
| dict thread | **51 %** (34 % `Router_Pop` on an empty inbox, 17 % `cpu_relax`) | 49 %: the insert's linear probing 34 %, the rest |
| producer | **26 %** (15 % `refill`'s `cpu_relax`, 7 % the ring's head, 4 % `free_pop_many`) | 74 % |

Both roles are idle a quarter to a half of the time, the wire is at 61 %, and the Router counts **10.5 M blocks
parked for a full destination and 1.34 M receives left unposted for want of a free block** (0 in every PROBE, 0
in the 1-node FILL).  A dict thread that waits half of the time inserts at 23 M/s when it works -- faster than
the 14.8 M/s of the 1-node FILL, where all 160 sit at the random-write wall together: the wall is not reached,
the points do not arrive.

Where the pool is: a node holds 62943 blocks of 262 KB; its producers cache 32 each (3040), its 160 inboxes hold
up to 64 each (10240), 32 are posted receives and 12 in flight, and the rest is sealed blocks parked because
their peer has its 4 credits in flight.  A send completes when the peer's matching receive does, and the peer
posts receives from its stash off the free ring; the senders leave a floor of `n_recv` = 32 blocks in that ring
and the service takes the last one.  The first reading was that the outbound side starves the inbound one: every
node's pool gone into blocks that wait for credit, no node with a block to receive with, no send completing.
**The intervention refutes it as stated**: a binary whose senders leave 8192 blocks in the ring (2.1 GB, 13 % of
the pool) runs FILL at exactly the same rate and still leaves receives unposted 2.6 M times.  A reserve that
large can only be eaten by the inbound path itself -- received blocks dispatched into inboxes and not yet given
back -- and the inboxes can hold 10240.  So the picture is **a burst regime**: the inflow to a node exceeds what
its dict threads take at moments (they are then at the random-write wall together), the inboxes fill, the ring
empties, receives go unposted, the peers' sends stall and their producers with them; then the inboxes drain and
the dict threads idle.  The averages -- dict threads idle 51 %, producers 26 %, wire at 61 % -- are what the
oscillation leaves.  How the knobs move it:

| FILL, 4 nodes, plateau a node | |
| --- | --- |
| 95/160 (225G) | 1.53 (1.52 on the repeat) |
| 63/192, fewer producers (225G) | **1.69** |
| 159/352, more producers (225G) | 0.60 |
| `--block 65536`, 4x the bytes a block (225G) | 0.98 |
| `--credit 16`, 4x the blocks in flight a peer (225G) | 1.29 |
| `--n-recv 128`, 4x the posted receives (225G) | 1.53 |
| 95/160 (900G) | 1.82 |
| `--inbox 256`, 4x the blocks a receiver may hold (900G) | 1.50 plateau, 1.25 exact against 1.82 / 1.27: nothing, and 2.1 M receives left unposted against 1.3 M |
| `--n-recv 128 --credit 8` (900G) | 1.81 plateau, 1.34 exact: nothing, 2.9 M receives left unposted |
| `--credit 2`, half the blocks in flight a peer (900G) | 1.79 plateau, 1.35 exact: nothing for FILL; PROBE loses 11 % (1.76 a node, the wire at 21.4 GB/s: 6 blocks in flight no longer cover the round trip) |
| **the senders' floor raised from 32 to 8192 blocks** (a one-line diagnostic build, `build-floor/`, not committed) | 1.81 plateau, 1.35 exact: **nothing**, and still 2.6 M receives left unposted |
| `--inbox 8`, an eighth of the blocks a receiver may hold (900G) | 1.80 plateau, 1.35 exact: nothing, 3.0 M receives left unposted; PROBE unchanged (24.1 GB/s) |

**Every flow-control knob is flat**: inboxes of 8, 64 or 256 blocks, credit 2, 4, 8 or 16 a peer, 32 or 128
posted receives, a 32- or 8192-block reserve -- FILL stays at 1.80 +- 0.02 a node at 900G, and the receives left
unposted stay in the millions whatever the reserve or the inbox, which says the blocks are neither only parked
outbound nor only sitting in inboxes: the pool is short at both ends by turns.  What FILL does respond to is the
team (63 producers 1.69, 95 producers 1.53, 159 producers 0.60, at 225G) and the block size (both 65 KB and 1 MB
blocks lose a third against 262 KB).  PROBE never enters the regime: its receivers take a block faster than the
wire fills one, and no PROBE ever left a receive unposted.  The 2-node FILL is the worst case, 1.07 a node.  **The
mechanism is open**; the burst picture above is the reading that fits, not a measurement, and the two things
still to look at are the service thread's turn in FILL (163 parked lists to retry, most of them full inboxes
whose head is a remote cache line, plus 62 K messages/s each way) and a per-destination flow control that
would keep a slow inbox from taking the node's inbound blocks.

## 5. Two nodes: 2.25 a node at 73 % of the wire

At 2 nodes half the points leave, so the wire would allow 24.3 / 8 = 3.0 G points/s a node, and the node alone
does 2.6; the run does 2.25 (4.16 exact), the wire at 18.3 GB/s.  One peer means 4 x 262 KB = 1 MB in flight
at most: at 18 GB/s that is 57 us, about a rendezvous's round trip, so the credit window and not the wire is the
2-node bound (at 4 nodes the three peers' windows add up).  It is not: `--credit 16` (4 MB in flight) gives 2.31 a node and 18.4 GB/s, `--block 65536` (16 MB) 2.31 and 17.3.
Whatever caps one peer pair at about 18.4 GB/s each way is below the Router: session 14's pure `router_bench` at
np 2 stopped at 85 % of `ucx_perftest`'s 23.5 GB/s on the same port, and this is 78 % of it with the dictionary
behind.  Open; the UCX side (`UCX_RNDV_SCHEME`, `UCX_TLS=dc_x`, one QP a peer) is where to look, and it stops
mattering from 4 nodes up, where the peers' shares add up to the port.

## What to change

1. **Fewer bytes a point on the wire.**  PROBE's ceiling on any grdix job is 24.3 GB/s / (16 B x (N - 1) / N);
   the destination's bits in the key are redundant once the block is addressed, and `val` is an n-bit preimage.
   Packing a point into 12 bytes (or 8, with the check bits cut to what the shard keeps) is the one change that
   raises the 4-node PROBE, up to the dict threads' 2.6 a node -- and it is what makes a 64-node run worth having.
2. **FILL's quarter is not a knob.**  Do not sweep `--inbox`, `--credit` or `--n-recv` for it again: 900G runs
   at every setting, and a build with a 256x larger reserve, all give 1.80 a node.  The next measurement is the
   service thread's own turn during FILL (its FILL profile was missed twice here: the profiler fired on PROBE),
   and the Router's idle-turn counter in the round report, which would say whether the service thread has slack
   in FILL as it has in PROBE; the next intervention is per-destination credit, so that a full inbox costs its own
   destination and not the node's inbound path.
3. **The service thread's core to itself** on the 512-thread team: the placer's SMT sibling costs the 4-node
   PROBE 5 %, and any faster wire or smaller message makes the service thread the bound.  Its own 35 % is
   scanning 163 parked lists a turn, most of them empty.
4. **Keep `--block` at 16384 points** (the auto choice here); 4096 halves the multi-node PROBE.
5. **Set `LD_LIBRARY_PATH` to the Guix Open MPI's lib** in every grdix run script, or give the demo a RUNPATH:
   the module does not, and the system 5.0.7 `libmpi.so.40` is what loads otherwise.

## Reproducing

```bash
oarsub -l {cluster='grdix'}/nodes=4,walltime=3:0:0 -q production -r ...          # 4 hosts, one rank each
ssh grdix-4; cd ~/mitm-grdix/build/session20
NP=4 N=42 P=95 R=160 RAM=900G TAG=np4_95_160 ./run20.sh                          # the reference, 10 min
python3 summ20.py np4_95_160; python3 live20.py 5 np4_95_160; python3 ib20.py np4_95_160 10
./batch20.sh                                                                     # the whole series, 90 min
# on grdix-5 while a run's PROBE (or FILL) is on: the per-role profiles
./waitprof20.sh TAG 160 95    # or waitprofF20.sh for FILL
```

Two lessons on the way.  **Never overwrite a running shell script** (bash reads it incrementally: an `scp` over
`run20.sh` while its `mpirun` ran made bash resume mid-word in the new file, run the demo a second time without
`mpirun`, and truncate the reference's log -- copy to a new name and `mv`).  And **`pkill -f` matches the shell
that runs it** when the pattern's text is on that shell's command line: `pkill -f '[i]bsample[.]sh'` in one ssh,
the `nohup bash ibsample.sh` in another.
