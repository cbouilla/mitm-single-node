# Parallel meet-in-the-middle

Given `f` and `g`, find a **claw** (`f(x0) == g(x1)`) or a **collision**
(`f(x0) == f(x1)`, `x0 != x1`), using van Oorschot–Wiener parallel collision search
(PCS) over a distributed dictionary of distinguished points.

The library is header-only C++17, in `include/`.  `examples/` holds one driver per
cipher: double-Speck64 (the current focus), double-DES, double-AES and SHA256.

There is **one engine**, and it is MPI + OpenMP.  One MPI rank per node; inside a
rank, thread 0 does all the MPI, some threads own a shard of the dictionary
("inserters") and the rest walk trails ("walkers").  A single rank with one walker
and one inserter is the degenerate sequential case — there is no separate sequential
engine to fall back on.

## Build

```bash
cmake -S . -B build && make -C build
```

Requires **MPI**, **OpenMP**, **OpenSSL** (headers; the DES example checks itself
against it) and **hwloc 2.x** (headers, `libhwloc-dev`): the threads are placed over
the NUMA nodes of the rank with it, and there is no fallback without it.  Binaries
land in `build/examples/`.

Two things to know before trusting any number that comes out:

- No `CMAKE_BUILD_TYPE` is set, so there is no `-O` flag and `NDEBUG` is off: the
  default build is unoptimized with live `assert()`s.  Add `-DCMAKE_BUILD_TYPE=Release`
  to benchmark.
- Everything compiles with `-march=native`, so binaries are **not portable across CPU
  generations**.  On a cluster, build on the same architecture as the compute nodes.
  The SIMD path (AVX-512 / AVX2 / scalar) is picked by `#ifdef` guards in
  `include/types.h` that `-march=native` triggers; the probes in `config/*.cmake` set
  `HAVE_*` variables that nothing currently reads.

## Run

Every driver takes the same options (`--help` lists them all).  `--ram` is
**mandatory** and accepts human units:

```bash
mpirun -np 4 --bind-to none build/examples/double_speck64_demo \
       --n 32 --ram 4G --inserters-per-node 2 --alpha 2.45 --beta 8
```

- `--n` problem size in bits (small == easy), `--seed` (0 == draw one and broadcast it)
- `--ram` dictionary bytes **per node**; `--alpha`, `--beta`, `--difficulty` tune
  the proportion of distinguished points and the length of a round
- `--nrounds` gives up after that many versions of the mixing function
- `--walkers-per-node` / `--inserters-per-node` set the thread layout; by default the
  walkers fill whatever affinity mask the launcher handed the rank, which is why
  `--bind-to none` matters.  The shards are pinned round-robin over the rank's NUMA
  nodes, so make `--inserters-per-node` a multiple of their count (rank 0 warns
  otherwise); `--no-bind` disables pinning
- the queue and buffer sizes (`--walker-queue`, `--buffer`, `--chunk`, ...) are the
  tuning knobs for the communication path

On a CEA/TGCC-style cluster the launcher is `ccc_mprun` under an `MSUB` batch script;
see `tgcc.sh`.

## Layout

```
include/
  types.h         SIMD vector types, selected by -march=native
  tools.hpp       PRNG (TRIVIUM), timing, human-readable numbers
  problem.hpp     the interface a cipher implements: f, g, is_good_pair, vfg
  parameters.hpp  Options (user knobs, all defaulted) and Parameters (derived once from them + RAM budget + problem size; data only)
  spsc.hpp        wait-free single-producer/single-consumer queue
  router.hpp      the Router, standalone
  comm.hpp        queues, bulk DP buffers, the counter enum, the control-channel payloads, ThreadContext and SharedContext
  walker.hpp      the walker thread; walking trails, turning a dictionary hit into a collision
  inserter.hpp    the inserter thread; PcsDict, the dictionary shard it builds and probes (held in SharedContext::shards)
  controller.hpp  rank 0's view: startup banner, rounds, pacing, statistics, when to stop
  engine.hpp      the comm thread (CommThread) and run(): preflight, the OpenMP team, the round loop
  mitm.hpp        umbrella: problem wrappers, claw_search(), collision_search()
  benchmark.hpp   f/g throughput per rank and across ranks, for the *_bench drivers
  naive/          the naive all-to-all MITM.  NOT PORTED, not built (see below)
examples/
  driver.hpp             command line + MPI startup, shared by every example
  <cipher>_problem.hpp   f, g, and a planted golden pair
  <cipher>_demo.cpp      run the attack
  <cipher>_bench.cpp     measure f evaluations per second
```

A new cipher is a `*_problem.hpp` implementing `f`, `g` and the domain against
`AbstractClawProblem` (or `AbstractCollisionProblem`), plus a ~20-line driver.
Vectorization is opt-in per problem: set `vlen` and provide `vfg` using the `v32` /
`v64` types from `types.h`.  With `vlen == 1` the engine calls `f`/`g` directly and no
`vfg` is needed.

Recurring notation: `n` = domain bits, `m` = range bits, `w` = dictionary slots,
`theta` = proportion of distinguished points, DP = distinguished point.

`include/naive/` and `examples/naive_double_speck64_demo.cpp` are a *different*
engine — the naive all-to-all MITM.  They still assume the old topology (several
ranks per node, split into senders and receivers) and do not compile; they are kept
as the starting point for porting that baseline back.

## Testing

CTest runs the Router's test suite only (`ctest --test-dir build`).  The
engine's demos self-check: each plants a golden pair, and both the problem wrappers and
the search entry points `assert` their way to it (`assert(pb.f(x0) == pb.g(x1))`).
Running `double_speck64_demo` on a small `--n` is the closest thing to a smoke test.
