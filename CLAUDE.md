# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A parallel meet-in-the-middle (MITM) attack framework for symmetric ciphers: given `f` and `g`, find **claws** (`f(x0) == g(x1)`) or **collisions**. The core is a header-only C++17 library in `include/`; `examples/` holds runnable drivers per cipher (double-Speck64 — the current focus, double-DES, double-AES, SHA256).

**There is exactly one engine, and it is MPI + OpenMP** (van Oorschot–Wiener parallel collision search over a distributed dictionary). The sequential and "vector sequential" engines were deleted; a single rank with one walker and one inserter is the degenerate sequential case. Do not reintroduce an engine-selection template parameter — `claw_search(pb, nbytes_memory, opts, prng)` and `collision_search(...)` take no engine.

## Build

```bash
cmake -S . -B build && make -C build
```

- C++17, compiled with `-march=native` — binaries are **not portable across CPU generations**, so on a cluster build on the same arch as the compute nodes. The SIMD path (AVX-512 / AVX2 / scalar) is selected by `#ifdef __AVX512F__ / __AVX2__` guards in `include/types.h` that `-march=native` triggers; the `config/*.cmake` probes set `HAVE_*` vars that are currently **unused**, and `config/neon.*` isn't wired into the root `CMakeLists.txt`.
- **No `CMAKE_BUILD_TYPE` / `-O` flag is set and `NDEBUG` is off**, so default builds are unoptimized with live `assert()`s. Keep this in mind before trusting benchmark numbers.
- Requires: **MPI**, **OpenMP**, **OpenSSL** (dev headers — the DES example self-checks against it).
- Binaries land in `build/examples/`. After changing CMake files, a clean rebuild may need `rm -rf build`.
- `examples/CMakeLists.txt` defines one `mitm_example(name [extra sources])` function; every driver links MPI + OpenMP + Threads the same way.

## Run

All drivers share one option parser (`examples/driver.hpp`, `--help` lists everything). **`--ram` is mandatory** (human units, e.g. `--ram 4G`): the search entry points take the byte budget per node as an explicit argument, and the *library* aborts on 0 — `driver.hpp` only parses the flag, it checks nothing. The bench drivers need no `--ram`.

```bash
mpirun -np 4 --bind-to none build/examples/double_speck64_demo \
    --n 32 --ram 4G --inserters-per-node 2 --alpha 2.45 --beta 8
```

Topology (the `Parameters` constructor in `include/parameters.hpp`): **one MPI rank per node**, `MPI_THREAD_FUNNELED`. Inside a rank, thread 0 is the comm thread (and, on rank 0, the controller), threads `1..inserters_per_node` own a dictionary shard each, the rest walk trails. Walkers default to filling the inherited affinity mask, hence `--bind-to none`. `--nrounds` bounds `max_versions` so a search can report "not found" instead of running forever. On a CEA/TGCC-style cluster the launcher is `ccc_mprun` under an `MSUB` batch script (`tgcc.sh`).

Quick smoke test (a few seconds):

```bash
mpirun -np 1 -bind-to none build/examples/double_speck64_demo --n 20 --ram 256K --difficulty 0.1 --walkers-per-node 6
```

## Layout & conventions

`include/` is flat — one engine, no `mpi/` or `sequential/` subdirectory:

- `types.h`, `tools.hpp`, `dict.hpp`, `problem.hpp` — SIMD types, PRNG/timing, `PcsDict`, the abstract problem interface.
- `counters.hpp` — per-round diagnostic tallies + HyperLogLog. **Round-scoped only**: threads accumulate, the comm thread merges and reduces them, and all-time totals plus every `printf` live in the controller.
- `parameters.hpp` — data only. `Options` is every user knob, all with a working default, no methods, no MPI. `Parameters : Options` is what the engine needs (topology, thread layout and placement as the `thread_cpu` table, dictionary size, difficulty), derived **once** by its constructor from `(opts, nbytes_memory, n, m)` inside `run_engine()` and const from then on. No method but the constructor. There is no `setup()`/`finalize()` protocol and no ordering to respect: `claw_search` takes `const Options &` and builds its own `Parameters`. Also `DP` and the MPI tags. The startup report is `Controller::banner()`; `benchmark()` (needs only `Options`, no RAM) lives in `benchmark.hpp`.
- `walker.hpp` / `inserter.hpp` — the two worker threads; `walker.hpp` also holds the trail code (`walk`, `walk_nolen1`, `resolve_collision`). `comm.hpp` / `spsc.hpp` — queues, bulk DP buffers, control channel, per-thread state. `controller.hpp` — rank 0's rounds and statistics. `engine.hpp` — `PcsNode` and `run_engine()`.
- `mitm.hpp` — umbrella: the problem wrappers and `claw_search` / `collision_search`.

**`PROTOCOL.md` is the bible of the engine.** It specifies every message exchanged between nodes (`TAG_POINTS`, `TAG_CONTROL`, the collectives, their wire layouts) and between threads (the SPSC queues, the collision queue, the per-thread `state` machine, the published tallies), plus every synchronization step of a round (start, steady state, the seven-step drain, epilogue, termination) and the loss/deadlock guarantees. **It must be kept in sync with the code**: any change to `comm.hpp`, `spsc.hpp`, `engine.hpp`, `controller.hpp`, `walker.hpp`, `inserter.hpp` or the tags/`DP` layout in `parameters.hpp` that touches a message, a queue, a state transition or the round lifecycle updates `PROTOCOL.md` in the same commit. If you find the code and the document disagreeing, say so rather than silently adjusting either.

`examples/driver.hpp` is *not* part of the library: it is the shared CLI + MPI startup for the drivers (`mitm::init`), and lives next to them.
- `naive/` — the naive all-to-all MITM. A **different** engine, **not ported and not built**: it still refers to `params.role`, `n_send`, `n_recv`, `recv_per_node`, `BCast_result`. Kept deliberately as the starting point for that work, along with `examples/naive_double_speck64_demo.cpp`. Don't try to "fix" it in passing.

Other conventions:

- A new cipher = a `*_problem.hpp` implementing `f`, `g`, domain against the abstract interface + a ~20-line driver.
- Vectorization is opt-in per problem: set `vlen` and provide `vfg`/`vf` using `v32`/`v64` from `types.h`. **With `vlen == 1` there is no `vfg` to provide** — the wrappers in `mitm.hpp` call `f`/`g` directly under `if constexpr`. The abstract problem classes deliberately *declare* `f`/`g`/`vfg` without defining them, so a missing member is a compile or link error rather than a function that silently returns 0 (which is what used to happen).
- Cipher primitives are C (`sha256.c`, `aes.c`, `usuba_des.cpp`); the framework is C++17.
- Notation: `n` = domain bits, `m` = range bits, `w` = dict slots, `theta` = distinguished-point proportion, DP = distinguished point.
- **Don't commit or hand-edit `build/`** — local, gitignored, out-of-source.

## Gotchas

- **`collision_search` is incomplete and always was.** `CollisionWrapper` in `mitm.hpp` (formerly `ConcreteCollisionProblem`, which carried a commented-out `assert(0); // not ready yet`) plateaus for some seeds and never retires the golden pair. Verified to fail identically on the deleted sequential engine, so it is the wrapper, not the engine — see the comment on the class. `claw_search` is the path that works.
- No CTest/unit-test suite. The demos self-check with `assert()` against a planted golden pair; running `double_speck64_demo` on a small `--n` is the smoke test.
- `theta == 1` ("zero difficulty") makes `start_chain` spin forever, because every point is distinguished and it refuses to start from one. `Controller::banner()` prints a loud warning; heed it.
- Percentages in the round report (`% probe failure` etc.) divide reduction totals by the DP count assembled from periodic reports. The two are collected differently, so a figure slightly over 100% is normal, not a bug.

## Workflow

Solo research repo (single remote, single author) — no PR process. Commit directly to feature branches (e.g. `omp_reboot`); keep it informal.
