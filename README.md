# Parallel meet-in-the-middle

Given `f` and `g`, find a **claw** (`f(x0) == g(x1)`) or a **collision**
(`f(x0) == f(x1)`, `x0 != x1`) over a distributed dictionary.

The library is header-only C++17, in `include/`.  `examples/` holds one driver per
cipher: double-Speck64 (the current focus), double-DES, double-AES and SHA256.

**One engine runs today: the direct, exhaustive meet-in-the-middle, built on the
Router.**  MPI + OpenMP, one rank per node.  Inside a rank the team is one service
thread doing all the MPI, `--dicts-per-node` threads owning a shard of the dictionary
each, and `--producers-per-node` threads evaluating the functions.  A round is two
phases: `f` fills the dictionary on a chunk of the domain, then `g` probes it on the
whole domain; `ceil(2^n / (fill * w))` rounds cover the domain, and finding nothing is
a proof of absence.  `PROTOCOL.md` §8 is its specification.

The **Router** (`include/router/`, spec `router.3`, `man -l router.3`) is the transport
underneath: it carries pairs of `u64` from sender threads to globally numbered receiver
threads across ranks, and it owns thread placement.

Parallel collision search (**PCS**, van Oorschot–Wiener, over a dictionary of
distinguished points) used to be the second engine.  It is **disconnected**: its files
and the older core they sit on are still in the tree, no target compiles them, and
`--engine pcs` is refused.  They are the record from which PCS will be rebuilt on the
Router.

## Build

```bash
cmake -S . -B build && make -C build
```

Requires **MPI**, **OpenMP**, **OpenSSL** (headers; the DES example checks itself
against it) and **hwloc 2.x** (headers, `libhwloc-dev`): the Router places the threads
with it, and there is no fallback without it.  Binaries land in `build/examples/`.

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
**mandatory** for the demos and accepts human units.  The exit status says whether the
golden pair was found.

```bash
mpirun -np 4 --bind-to none build/examples/double_speck64_demo \
       --n 32 --ram 4G --dicts-per-node 2 --producers-per-node 30
```

- `--n` problem size in bits (small == easy), `--seed` (0 == draw one and broadcast it)
- `--ram` dictionary bytes **per node**; `--fill` the fill ratio, so `fill * w` entries
  per round; `--nrounds` gives up after that many rounds
- `--producers-per-node` / `--dicts-per-node` size the team.  Producers default to
  filling whatever affinity mask the launcher handed the rank, which is why
  `--bind-to none` matters
- **the Router pins the threads and forms its groups itself**, and prints both the
  planned and the measured layout.  `--no-bind` turns pinning off, which is what two
  ranks sharing a host need (the Router refuses to pin overlapping masks);
  `--cache-level` and `--group` are its group knobs
- `--block`, `--swc`, `--n-recv`, `--inbox`, `--sweep`, `--credit` tune the transport;
  `router.3` documents each of them

A smoke test that runs in well under a second:

```bash
mpirun -np 1 --bind-to none build/examples/double_speck64_demo --n 20 --ram 256K \
       --producers-per-node 4 --dicts-per-node 3
```

On a CEA/TGCC-style cluster the launcher is `ccc_mprun` under an `MSUB` batch script;
see `tgcc.sh`.

## Layout

```
include/
  types.h            SIMD vector types, selected by -march=native
  tools.hpp          PRNG (TRIVIUM), timing, hashing, human-readable numbers, Point
  problem.hpp        the interface a cipher implements: f, g, is_good_pair, vfg
  parameters.hpp     Options: every user knob, all defaulted, the Router's among them
  benchmark.hpp      f/g throughput per rank and across ranks, for the *_bench drivers
  direct_common.hpp  the direct engine: parameters, counters, the shared tallies
  direct_dict.hpp      the dictionary shard and the thread that fills and probes it
  direct_producer.hpp  the thread that evaluates the phase's function and pushes
  direct.hpp           the problem wrappers, the service round, the epilogue, run()
  router/            the Router: ring, common, placement, connect, workers, service
                     (router.hpp is the umbrella; router.3 at the root is its man page)
  pcs_common.hpp, walker.hpp, inserter.hpp, pcs.hpp    PCS -- disconnected
  comm.hpp, spsc.hpp, controller.hpp, engine.hpp, placement.hpp
                     the older core PCS sits on -- disconnected, does not compile
examples/
  driver.hpp             command line + MPI startup, shared by every example
  router_driver.hpp      the same for the Router's own two programs
  <cipher>_problem.hpp   f, g, and a planted golden pair
  <cipher>_demo.cpp      run the attack
  <cipher>_bench.cpp     measure f evaluations per second
  router_test.cpp        the Router's test suite
  router_bench.cpp       the Router's throughput benchmark
```

A new cipher is a `*_problem.hpp` implementing `f`, `g` and the domain against
`AbstractClawProblem` (or `AbstractCollisionProblem`), plus a ~20-line driver.
Vectorization is opt-in per problem: set `vlen` and provide `vfg` using the `v32` /
`v64` types from `types.h`.  With `vlen == 1` the engine calls `f`/`g` directly and no
`vfg` is needed.

Recurring notation: `n` = domain bits, `m` = range bits, `w` = dictionary slots,
`theta` = proportion of distinguished points, DP = distinguished point.

## Documents

- `router.3` — the Router's interface, the authoritative reference for it.
- `PROTOCOL.md` — §8 specifies the direct engine; §1–§7 are the older core and PCS,
  kept as the record to rebuild PCS from.
- `PROBLEM.md` — the measurements that showed the single comm thread was the
  bottleneck, which is why the Router exists.
- `FINDINGS.md` — one section per profiling session on Grid'5000.

## Testing

```bash
ctest --test-dir build
```

Five tests: the Router's suite at one, two and four ranks, and the direct engine at one
and two ranks.  The engine's tests are `double_speck64_demo` on a small `--n` against a
planted golden pair — the demos `assert` their way to it
(`assert(pb.f(x0) == pb.g(x1))`) and their exit status reports whether they found it.
