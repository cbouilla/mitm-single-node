# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A parallel meet-in-the-middle (MITM) attack framework for symmetric ciphers: given `f` and `g`, find **claws** (`f(x0) == g(x1)`) or **collisions**. The core is a header-only C++17 library in `include/`; `examples/` holds runnable drivers per cipher (double-Speck64 — the current focus, double-DES, double-AES, SHA256).

**There is exactly one engine, and it is MPI + OpenMP** (van Oorschot–Wiener parallel collision search over a distributed dictionary). The sequential and "vector sequential" engines were deleted; a single rank with one producer and one dict thread is the degenerate sequential case. Do not reintroduce an engine-selection template parameter — `claw_search(pb, nbytes_memory, opts, prng)` and `collision_search(...)` take no engine.

## Build

```bash
cmake -S . -B build && make -C build
```

- C++17, compiled with `-march=native` — binaries are **not portable across CPU generations**, so on a cluster build on the same arch as the compute nodes. The SIMD path (AVX-512 / AVX2 / scalar) is selected by `#ifdef __AVX512F__ / __AVX2__` guards in `include/types.h` that `-march=native` triggers; the `config/*.cmake` probes set `HAVE_*` vars that are currently **unused**, and `config/neon.*` isn't wired into the root `CMakeLists.txt`.
- **No `CMAKE_BUILD_TYPE` / `-O` flag is set and `NDEBUG` is off**, so default builds are unoptimized with live `assert()`s. Keep this in mind before trusting benchmark numbers.
- Only the MPI **C** API is used, so `find_package(MPI REQUIRED COMPONENTS C)` and the drivers link `MPI::MPI_C`.
  Don't ask for the `CXX` component: it is not the C++ bindings (that is `MPICXX`), it is FindMPI probing the
  `mpicxx` wrapper, which some clusters ship broken or not at all — a spurious `missing: MPI_CXX_FOUND`.
- Requires: **MPI**, **OpenMP**, **OpenSSL** (dev headers — the DES example self-checks against it), **hwloc 2.x** (dev headers, `libhwloc-dev`): it places the threads over the NUMA nodes and is a hard requirement with no fallback. If CMake can't find it, `module load hwloc` or `HWLOC_ROOT=/path cmake ...`.
- Binaries land in `build/examples/`. After changing CMake files, a clean rebuild may need `rm -rf build`.
- `examples/CMakeLists.txt` defines one `mitm_example(name [extra sources])` function; every driver links MPI + OpenMP + Threads + hwloc the same way.

## Run

All drivers share one option parser (`examples/driver.hpp`, `--help` lists everything). **`--ram` is mandatory** (human units, e.g. `--ram 4G`): the search entry points take the byte budget per node as an explicit argument, and the *library* aborts on 0 — `driver.hpp` only parses the flag, it checks nothing. The bench drivers need no `--ram`.

```bash
mpirun -np 4 --bind-to none build/examples/double_speck64_demo \
    --n 32 --ram 4G --dicts-per-node 2 --alpha 2.45 --beta 8
```

Topology (`include/placement.hpp`): **one MPI rank per node**, `MPI_THREAD_FUNNELED`. Inside a rank, thread 0 is the comm thread (and, on rank 0, the controller), threads `1..dicts_per_node` own a dictionary shard each, the rest walk trails. Producers default to filling the inherited affinity mask, hence `--bind-to none`. **The engine assumes one MPI rank per NUMA node** and rank 0 warns loudly otherwise. Placement is cache- and SMT-aware (hwloc): the mask's **cores** are cut into **one thread group per dict thread**, each inside one cache domain — the lowest level shared by several *cores*, L3 on both Xeon and EPYC, detected rather than hardcoded and overridable with `--cache-level`. Groups own whole cores so that a core's SMT siblings stay in one group. Every thread then takes the **emptiest core of its group**, service threads first, which gives the comm thread and each dict thread a core of its own and hands the siblings they leave to producers — a producer fills the issue slots a thread waiting on queues or DRAM leaves idle, and two waiting threads must never share a core (`PROBLEM.md` §3: ~10% to them, +66% to the producers). **A producer's group IS the dict thread whose collisions it resolves**, so `--dicts-per-node` sets the group count: make it at least the number of cache domains (rank 0 warns otherwise) and never more than `--producers-per-node` (refused: every collision queue needs a producer). `--no-bind` disables pinning and falls back to `walker_index % dicts_per_node`. Each thread records the CPU and NUMA node the kernel reports after pinning in its `ThreadContext` (`cpu`, `numa_node`); a thread that isn't where it should be aborts the run before anything is allocated. Other topologies are testable on a laptop with `HWLOC_SYNTHETIC="node:1 l3:4 core:4 pu:1"` set on the *application* (`mpirun ... env HWLOC_SYNTHETIC=... build/examples/...`), not on `mpirun`, whose own mapper would then see the fake machine. `--nrounds` bounds `max_versions` so a search can report "not found" instead of running forever. On a CEA/TGCC-style cluster the launcher is `ccc_mprun` under an `MSUB` batch script (`tgcc.sh`).

Quick smoke test (a few seconds):

```bash
mpirun -np 1 -bind-to none build/examples/double_speck64_demo --n 20 --ram 256K --difficulty 0.1 --producers-per-node 6
```

Add `--dp-len-bits 1` to that command for the resolver's torture test: it ships every trail length as
saturated, and a length unknown on the wire is stored unknown too, so **every** candidate takes the
"both lengths unknown" path (`PROTOCOL.md` §3.2) instead of the ~0% it sees at the default width. It
must still find the golden pair, with `walk-noncolliding` no worse than the default run. It is a
correctness test only, never a measurement: with every length unknown the dictionary's "keep the trail
at least as long" rule always overwrites, so the length bias of the collisions found is gone and
`COLLIDING_LEN_MIN/MAX` shift. `sha2_claw_demo` is the only driver with `vlen == 1`, hence the only one
that exercises the scalar resolver rather than `VecResolver`.

## Layout & conventions

`include/` is flat — one engine, no `mpi/` or `sequential/` subdirectory:

- `types.h`, `tools.hpp`, `problem.hpp` — SIMD types, PRNG/timing/hashing/formatting, the abstract problem interface.
- `placement.hpp` — **the only file that knows about CPUs**: hwloc (NUMA nodes, cache domains, cores), the partition of the affinity mask into one thread group per dict thread, `pin_to_cpu`, and the two reports (the plan, printed by the banner; the measured layout, printed once the team is pinned). `Parameters` holds one `Placement`; nothing else in the library queries the topology.
- `parameters.hpp` — data only. `Options` is every user knob, all with a working default, no methods, no MPI. `Parameters : Options` is what the engine needs (topology, the `Placement` (thread layout, see above), dictionary size, difficulty), derived **once** from `(opts, nbytes_memory, n, m)` at the top of `run()` and const from then on. No method but the constructor. There is no `setup()`/`finalize()` protocol and no ordering to respect: `claw_search` takes `const Options &` and hands them to `run()`, which derives the `Parameters`. Also `DP` and the MPI tags. The startup report is the static `Controller::banner()`, printed by `run()` before the team exists; `benchmark()` (needs only `Options`, no RAM) lives in `benchmark.hpp`.
- `walker.hpp` / `inserter.hpp` — the two worker threads; `walker.hpp` also holds the trail code (`walk`, `measure_trail`, `walk_recorded`, `resolve_collision`), `inserter.hpp` also holds `PcsDict`, the dictionary shard a dict thread builds and probes, held in `SharedContext::shards` (`comm.hpp` only forward-declares it). `comm.hpp` / `spsc.hpp` — queues, bulk DP buffers, `enum counter`, the control-channel payload layouts (the channel's receives themselves live in `CommThread` and `Controller`), and the two context structs: `ThreadContext` (what one thread owns and the comm thread also touches: `role`, the measured `cpu` and `numa_node`, its `group`, `state`, the `u64 ctr[N_COUNTERS]` tallies — plain, no atomics — a worker's SPSC queue, and a producer's own HyperLogLog, also plain) and `SharedContext` (what all threads share: the round header, the `ctx`, `shards` and `coll_q` tables, the golden pair). The counter array is the wire layout of a progress report (deltas) and of the end-of-round `MPI_Reduce` (totals); it is round-scoped, and all-time totals plus every `printf` live in the controller. `controller.hpp` — rank 0's rounds and statistics. `engine.hpp` — `CommThread` (thread 0, all of it: routing, control channel, the round's state machine, epilogue, termination) and `run()`, the engine's one entry point: preflight, the OpenMP team and its round loop; what the threads share (`Parameters`, the `SharedContext`) and the `CommThread` are its locals.
- `mitm.hpp` — umbrella: the problem wrappers and `claw_search` / `collision_search`.

**`PROTOCOL.md` is the bible of the engine.** It specifies every message exchanged between nodes (`TAG_POINTS`, `TAG_END_ROUND`, `TAG_REPORT`, `TAG_SOLUTION`, the collectives, their wire layouts) and between threads (the SPSC queues, the per-dict thread collision queues, the per-thread `state` machine, the per-thread tallies), plus every synchronization step of a round (start, steady state, the seven-step drain, epilogue, termination) and the loss/deadlock guarantees. **It must be kept in sync with the code**: any change to `comm.hpp`, `spsc.hpp`, `engine.hpp`, `controller.hpp`, `walker.hpp`, `inserter.hpp` or the tags/`DP` layout in `parameters.hpp` that touches a message, a queue, a state transition or the round lifecycle updates `PROTOCOL.md` in the same commit. If you find the code and the document disagreeing, say so rather than silently adjusting either.

`examples/driver.hpp` is *not* part of the library: it is the shared CLI + MPI startup for the drivers (`mitm::init`), and lives next to them.
- `naive/` — the naive all-to-all MITM. A **different** engine, **not ported and not built**: it still refers to `params.role`, `n_send`, `n_recv`, `recv_per_node`, `BCast_result`. Kept deliberately as the starting point for that work, along with `examples/naive_double_speck64_demo.cpp`. Don't try to "fix" it in passing.

Other conventions:

- A new cipher = a `*_problem.hpp` implementing `f`, `g`, domain against the abstract interface + a ~20-line driver.
- Vectorization is opt-in per problem: set `vlen` and provide `vfg`/`vf` using `v32`/`v64` from `types.h`. **With `vlen == 1` there is no `vfg` to provide** — the wrappers in `mitm.hpp` call `f`/`g` directly under `if constexpr`. The abstract problem classes deliberately *declare* `f`/`g`/`vfg` without defining them, so a missing member is a compile or link error rather than a function that silently returns 0 (which is what used to happen).
- Cipher primitives are C (`sha256.c`, `aes.c`, `usuba_des.cpp`); the framework is C++17.
- Notation: `n` = domain bits, `m` = range bits, `w` = dict slots, `theta` = distinguished-point proportion, DP = distinguished point.
- **Where state lives**, per rank: private to a worker thread and never touched by the comm thread → a local of that thread's function; specific to one thread but touched by the comm thread (`role`, `state`, `ctr`, its SPSC queue) → its `ThreadContext`; common to all threads → the `SharedContext`; private to the comm thread (MPI buffers, requests, the controller) → a member of `CommThread`. Two refinements: **every counter, the comm thread's own included, goes in `ctx[tid]->ctr[]`, no exceptions**, so the node's tallies are one loop with no special case; and **the dictionary shards and the collision queues live in `SharedContext::shards` / `::coll_q`** although only their dict thread produces into them, because zeroing a shard may one day be collective and because a queue's consumers are a whole thread group. A producer's HyperLogLog follows the counter rule: one plain array per producer in its `ThreadContext`, merged by the comm thread at the end of a round.
- **Per-thread objects are built by their owning thread**, in `run()` right after `pin_to_cpu`: its `ThreadContext` into `shared.ctx[tid]` (whose constructor builds a worker's SPSC queue) and, for a dict thread, its `PcsDict` shard into `shared.shards[tid-1]`. The `CommThread` is the exception: `run()` builds it before the team, on the main thread, which becomes thread 0 (its MPI buffers are tens of kB, placement is irrelevant). Before the team exists, `run()` only derives the `Parameters` and builds the `SharedContext`, whose two tables are sized but empty. This is NUMA first touch — the shard's zero-fill is what places its pages — so don't move construction back into a constructor. Deliberately left for later: a NUMA-affine parallel flush of the shards between rounds (which threads share a shard's NUMA node is now known from `ThreadContext::numa_node`; the collective flush itself isn't written), and a generation bit in dictionary entries to skip the flush altogether (changes the "empty at round start" guarantee in `PROTOCOL.md`).
- **Don't commit or hand-edit `build/`** — local, gitignored, out-of-source.

## Gotchas

- **`collision_search` is incomplete and always was.** `CollisionWrapper` in `mitm.hpp` (formerly `ConcreteCollisionProblem`, which carried a commented-out `assert(0); // not ready yet`) plateaus for some seeds and never retires the golden pair. Verified to fail identically on the deleted sequential engine, so it is the wrapper, not the engine — see the comment on the class. `claw_search` is the path that works.
- No CTest/unit-test suite. The demos self-check with `assert()` against a planted golden pair; running `double_speck64_demo` on a small `--n` is the smoke test.
- `theta == 1` ("zero difficulty") makes `start_chain` spin forever, because every point is distinguished and it refuses to start from one. `Controller::banner()` prints a loud warning; heed it.
- The live one-line display and the decision to close a round rest on the progress reports, which are approximate (one-way, and read from the workers' tallies without synchronisation). The round report printed at the end of a round comes from the exact `MPI_Reduce`. Don't expect the two to agree to the unit.

## Workflow

Solo research repo (single remote, single author) — no PR process. Commit directly to feature branches (e.g. `omp_reboot`); keep it informal.
