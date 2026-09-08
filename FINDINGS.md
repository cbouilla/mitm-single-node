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
