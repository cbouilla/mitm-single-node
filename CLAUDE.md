# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A parallel meet-in-the-middle (MITM) attack framework for symmetric ciphers: given `f` and `g`, find **claws** (`f(x0) == g(x1)`) or **collisions**. Header-only C++17 in `include/`; `examples/` holds runnable drivers per cipher (double-Speck64 — the current focus, double-DES, double-AES, SHA256).

**One engine runs today: the direct engine, on the Router.** `direct_common.hpp`, `direct_dict.hpp`, `direct_producer.hpp`, `direct.hpp` (`PROTOCOL.md` §8), entry points `mitm::direct::claw_search(pb, nbytes_memory, opts, prng)` and `mitm::direct::collision_search(...)`. It is the exhaustive meet-in-the-middle: `ceil(2^n / (fill * w))` rounds of two phases, `f` fills the distributed dictionary on a chunk of the domain, then `g` probes it on the whole domain, matches resolved by the dict threads themselves, nothing dropped. **It is built on the Router alone** (`include/router/`, spec `router.3`) and uses none of the old core: no comm thread, no control channel, no controller, no SPSC queue, no scheme template. Per rank, one OpenMP team of `1 + dicts_per_node + producers_per_node` threads: thread 0 is the Router's service thread, the next ones are receivers owning a dictionary shard each, the rest are senders evaluating the phase's function. One direct phase is one Router round, so `Router_Reset` is the FILL→PROBE barrier the algorithm needs; the epilogue is one `MPI_Allgather` of the node's tallies, its `Router_Stats` and its golden pair, from which every node reads the same verdict. **The Router owns placement in autopilot mode**: the engine passes `ROUTER_GROUP_AUTO`, never pins a thread and never queries hwloc.

**PCS is disconnected, on purpose.** `pcs_common.hpp`, `walker.hpp`, `inserter.hpp`, `pcs.hpp` and the old core they sit on — `engine.hpp`, `comm.hpp`, `controller.hpp`, `spsc.hpp`, `placement.hpp` — are still in the tree but **no target compiles them and they are expected not to build**; `--engine pcs` is refused with a message. They are the record to rebuild PCS on the Router from (`PROTOCOL.md` §1–§7). Do not repair them in passing, and do not add them back to a target: rebuilding PCS on the Router is its own job. The sequential and "vector sequential" engines were deleted long before; a single rank with one producer and one dict thread is the degenerate sequential case.

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

All drivers share one option parser (`examples/driver.hpp`, `--help` lists everything). **`--ram` is mandatory** for the demos (human units, e.g. `--ram 4G`): the search entry points take the byte budget per node as an explicit argument, and the *library* aborts on 0 — `driver.hpp` only parses the flag, it checks nothing. The bench drivers need no `--ram`. `--engine` accepts `direct` only (the default); `--fill` is the dictionary fill ratio (default 0.5). Every driver's exit status says whether the golden pair was found.

```bash
mpirun -np 4 --bind-to none build/examples/double_speck64_demo \
    --n 32 --ram 4G --dicts-per-node 2 --producers-per-node 30
```

Topology: **one MPI rank per node**, `MPI_THREAD_FUNNELED`, and **the Router decides where every thread runs** — it pins them, forms its groups and prints both the planned and the measured layout itself. `--dicts-per-node` and `--producers-per-node` (0 == fill the inherited affinity mask) size the team; producers default to filling the mask, hence `--bind-to none`. `--no-bind` hands `opts.router.pin = false` to the Router, which is what two ranks sharing a host need (it refuses to pin overlapping masks). `--cache-level` and `--group` are the Router's group knobs, `--block`, `--swc`, `--n-recv`, `--inbox`, `--sweep`, `--credit` its transport knobs; `router.3` documents all of them, and the same flags exist in `examples/router_driver.hpp`. Other topologies are testable on a laptop with `HWLOC_SYNTHETIC="node:1 l3:4 core:4 pu:1"` set on the *application* (`mpirun ... env HWLOC_SYNTHETIC=... build/examples/...`), not on `mpirun`, whose own mapper would then see the fake machine. `--nrounds` bounds the round count so a search can report "not found" instead of running to exhaustion. On a CEA/TGCC-style cluster the launcher is `ccc_mprun` under an `MSUB` batch script (`tgcc.sh`).

Quick smoke test, 66 rounds of two phases over 2^20 preimages (under a second; `--nrounds 3` makes it give up instead, and exit nonzero):

```bash
mpirun -np 1 --bind-to none build/examples/double_speck64_demo --n 20 --ram 256K \
    --producers-per-node 4 --dicts-per-node 3
```

Two ranks on the one laptop, unpinned (8 spinning threads in all):

```bash
mpirun -np 2 --bind-to none build/examples/double_speck64_demo --n 20 --ram 256K \
    --producers-per-node 2 --dicts-per-node 1 --no-bind
```

`sha2_claw_demo` is the only driver with `vlen == 1`, hence the only one that exercises the scalar `veval` rather than the problem's `vfg`; `sha2_collision_demo` is the only one that exercises `CollisionWrapper` and its two-order `is_good_pair`.

## Layout & conventions

`include/` is flat, and it holds two things: the direct engine plus the Router (what builds), and the old core plus PCS (what does not, pending PCS's rebuild).

- `types.h`, `tools.hpp`, `problem.hpp` — SIMD types; PRNG/timing/hashing/formatting and `Point` (the two-word wire unit the Router carries: routing key + payload); the abstract problem interface.
- `parameters.hpp` — **`Options` and nothing else**: every user knob, all with a working default, no methods, no MPI, plus one `Router_Opts router` member holding the Router's own. A driver fills it and hands it to an entry point, which derives its parameters once. `benchmark.hpp` is `benchmark()` alone, the problem's f/g rate: `Options` and a problem, no RAM, no engine.
- `direct_common.hpp` / `direct_dict.hpp` / `direct_producer.hpp` / `direct.hpp` — **the direct engine**, `namespace mitm::direct`: its data (`direct::Params` with the team's shape, the domain, `per_round`, `n_rounds`, `check_bits`; the 7 counters and the round record's layout; `Tally` and `Shared`, the tallies plus the golden pair and the verdict); the linear-probing `DirectDict` and the dict thread that inserts in `FILL` and probes-and-resolves in `PROBE`; the producer that enumerates its piece of the domain and pushes; then the two wrappers (`ClawWrapper`, `CollisionWrapper`: `fill()`, `veval()`, `good()`, `self_test()`), the service round, the epilogue, all the printing, `run()` and `direct::claw_search` / `direct::collision_search`. `PROTOCOL.md` §8.
- `router/topology.hpp` — **the hwloc primitives**: `pin_to_cpu`, `lowest_shared_cache_level`, `topology_of_cpus`, `cores_by_domain`, `emptiest_core`, in `namespace mitm`. The one file that touches hwloc; the Router's `router/placement.hpp` and the old core's `placement.hpp` both build on it. The direct engine calls none of it — placement is the Router's.

The disconnected half, kept for PCS's rebuild and **expected not to compile**:

- `pcs_common.hpp` — the PCS scheme's data, in `namespace mitm::pcs`: `pcs::Params`, `Header` (`i`, `root_seed`), `enum counter`, `CollisionCandidate` / `CollisionQueue`, the HyperLogLog, `ThreadStats` / `RoundStats`, `Shared`, and the declaration of `pcs::Scheme`. `walker.hpp` / `inserter.hpp` — PCS's two worker threads; `walker.hpp` also holds the trail code (`walk`, `measure_trail`, `walk_recorded`, `resolve_collision`), `inserter.hpp` also holds `PcsDict`. `pcs.hpp` — the three problem wrappers, the definitions of `pcs::Scheme`'s functions and `pcs::claw_search` / `pcs::collision_search`.
- `comm.hpp` / `spsc.hpp` — the old core's `Parameters` (which moved here out of `parameters.hpp`, with `enum tags` and `POINT_WORDS`), the SPSC ring of `Point`s, the bulk point buffers, the control-channel payload layouts, and the two context structs `ThreadContext` / `SharedContext`. `controller.hpp` — rank 0's round mechanics and `layout_banner()`. `engine.hpp` — `CommThread<Scheme>` and `run<Scheme>()`. `placement.hpp` — the old core's thread plan, one group per dict thread.

**`PROTOCOL.md` is the bible of the engine.** §8 specifies the direct engine: the team and its Router roles, the phase sequence and why `Router_Reset` is the barrier the algorithm needs, the wire point and its two hashings, the epilogue's `MPI_Allgather` record, where state lives, the tallies, the loss guarantee. §1–§7 are the old core and PCS, kept as the record to rebuild PCS from. **§8 must be kept in sync with the code**: any change to the four `direct_*.hpp` files that touches the phase lifecycle, the point format, the routing, the record's layout, a barrier or the loss semantics updates it in the same commit. If you find the code and the document disagreeing, say so rather than silently adjusting either. `router.3` is the Router's own spec and is authoritative for everything the engine calls (`man -l router.3`).

`examples/driver.hpp` is *not* part of the library: it is the shared CLI + MPI startup for the drivers (`mitm::init`), and lives next to them.

**The Router** (`include/router/`, a directory of one header per concern; the umbrella `router.hpp`'s header comment is its summary and `router.3` at the root its man page, `man -l router.3`) is **the transport the direct engine runs on**: an MPI-style procedural interface around one RAII handle, the `Router_thread` that `Router_Init(role, group, comm, tag, lossy, opts)` returns to each thread by value (collective over every thread of every node, the same arguments everywhere **but `group`**, `ROUTER_DEFAULT_OPTS` for the defaults; the caller keeps it named for the team's whole life). **The Router owns placement**: by default it pins the threads and forms G groups, a group being senders and receivers in one cache domain, at most `opts.group_size` cores each — the receivers then the senders are round-robined over the domains into the least-filled group (`Router_group` / `Router_num_groups` / `Router_group_num_receivers` / `Router_group_receiver(i)` expose it, a producer pairs with the dict thread of its group; `Router_domain` / `Router_num_domains` / `Router_cpu` / `Router_numa_node` report where the kernel put a thread). `opts.pin = false` leaves the CPUs to the caller and groups by worker index; a per-thread `group` color takes the caller's groups à la `MPI_Comm_split`. When pinning it **refuses** a mask too small for the team or one shared with a co-hosted rank, so two ranks on one host need `--no-bind`. Every round call takes the handle: `Router_Push`, `Router_Close`, `Router_Grab`, `Router_Release`, `Router_Pop`, `Router_Test_drained` a worker's own (no `omp_get_thread_num()` on the hot path), `Router_Progress`, `Router_Test_quiescent`, `Router_Stats` the service thread's; `Router_Reset` is collective over the team (its two team barriers and the `MPI_Barrier` are inside it, so the caller writes none). The node's shared state (`Router_node`) is internal: thread 0's handle owns it and deletes it when it dies at the end of the parallel region, so everything is gone before `MPI_Finalize` with no scoping needed. It routes points from sender threads to globally numbered receiver threads across ranks: senders stage into shared per-destination blocks through a private write-combining line (PROBLEM.md §9), with **non-temporal stores** (`router_stream_copy`), so a block is never fetched into the writer's cache — measured 2.1x on a 256-core EPYC, where a plain `memcpy` pulled every staged line in from DRAM, half of it across the sockets; a block is a message: the service thread (thread 0, all MPI) sends it from block memory, receives into block memory and hands blocks to the receivers by id (the pool is first-touched by every thread of the team, a slice of blocks each, so its pages spread over the rank's NUMA nodes instead of all landing on the service's), copying nothing at all: a receiver reads a block in place between `Router_Grab` and `Router_Release`, `Router_Pop` being one point at a time on top of them; in lossy mode no call ever waits. Free blocks live in a bounded MPMC ring of block ids (Vyukov): a receiver pushes by a `fetch_add` ticket and never fails, a sender or the service pops a batch by one CAS and fails cleanly when empty. It shares only `tools.hpp` (for `Point`) and `topology.hpp` with the rest of the tree; its `RouterRing` is a minimal, non-templated copy of `spsc.hpp`, on purpose. The parts: `router.hpp` (umbrella + spec), `ring.hpp`, `common.hpp` (constants, `Router_Opts`, `Router_thread`, `Router_node`), `placement.hpp` (`RouterPlacement`), `connect.hpp` (`Router_Init`, `connect`, the queries), `workers.hpp` (the free ring, the sender and receiver paths), `service.hpp` (the service thread). `examples/router_test.cpp` is its test suite (`ctest --test-dir build`, or `mpirun -np 2 --bind-to none build/examples/router_test --senders 3 --receivers 2 --no-bind`; `--test NAME` for one, `groups`/`colors` cover the placement), `examples/router_bench.cpp` its benchmark (`--help`); both take their options from `examples/router_driver.hpp`. The laptop has 8 physical cores: keep its runs to 8 spinning threads in all (np 1 pins: `--senders 4 --receivers 3`; np 2 `--senders 2 --receivers 1 --no-bind`; np 4 needs 12, correctness only).

Other conventions:

- A new cipher = a `*_problem.hpp` implementing `f`, `g`, domain against the abstract interface + a ~20-line driver.
- Vectorization is opt-in per problem: set `vlen` and provide `vfg`/`vf` using `v32`/`v64` from `types.h`. **With `vlen == 1` there is no `vfg` to provide** — the wrappers in `direct.hpp` call `f`/`g` directly under `if constexpr`. The abstract problem classes deliberately *declare* `f`/`g`/`vfg` without defining them, so a missing member is a compile or link error rather than a function that silently returns 0 (which is what used to happen).
- Cipher primitives are C (`sha256.c`, `aes.c`, `usuba_des.cpp`); the framework is C++17.
- Notation: `n` = domain bits, `m` = range bits, `w` = dict slots, `theta` = distinguished-point proportion, DP = distinguished point.
- **Where state lives**, per rank (the direct engine): private to a thread and touched by nobody else — its Router handle, a dict thread's shard, a producer's buffers, thread 0's stats and record arrays → a local of that thread's scope; shared by the team — the tallies, the golden pair, the verdict → the one `Shared`, built before the team. `Shared::tally[tid].ctr[]` is one `u64[N_COUNTERS]` per thread on a cache line of its own, **plain, never atomic**: its owner writes it and thread 0 reads it after `Router_Reset`'s first team barrier, which is what makes the writes visible. The golden pair is the one exception (a mutex and an atomic flag), because any dict thread may find one at any moment.
- **Per-thread objects are built by their owning thread**, right after `Router_Init` — i.e. once the Router has pinned it. A dict thread's shard is the case that matters: its zero-fill is the NUMA first touch of every page, so don't move it before the team or into a constructor called elsewhere. Deliberately left for later: a NUMA-affine parallel flush of the shards between rounds, and a generation bit in dictionary entries to skip the flush altogether (which would change the "empty at the start of a FILL phase" guarantee in `PROTOCOL.md`).
- **Don't commit or hand-edit `build/`** — local, gitignored, out-of-source.

## Gotchas

- **PCS does not build, and that is the current state, not a bug.** `--engine pcs` is refused; `pcs.hpp` and the old core are unreferenced by every target. Don't fix them in passing.
- **`pcs::collision_search` was incomplete and always was**, and that is waiting for PCS's rebuild: its `CollisionWrapper` plateaus for some seeds and never retires the golden pair. `direct::collision_search` is exhaustive and works: it is the reference to debug it against when the time comes.
- CTest runs the Router's tests and the direct engine's (`router_np1/2/4`, `direct_np1/2`). The engine has no unit tests beyond that: the demos plant a golden pair and `assert()` their way to it, and their **exit status** is whether they found it — which is what the two `direct_*` tests check. With Open MPI, CMake defaults `MPIEXEC_PREFLAGS` to `--bind-to none`, or a rank's whole team lands on one core.
- **The live one-line display is rank 0's own node**, read from its tallies without synchronisation and scaled by the node count: approximate on purpose. The round report comes from the epilogue's `MPI_Allgather` and is exact. Don't expect the two to agree to the unit.
- The direct engine has **no early exit inside a phase**: a golden pair found during a phase stops the search at the end of that phase, not at once. Cutting a phase short would need a control channel of its own and can save at most one phase.
- A dict thread reads each delivered block **in place**, between `Router_Grab` and `Router_Release`. A block may hold fewer than `block_points` points, and a receiver that returns while holding one both prevents quiescence and aborts in `Router_Reset`.
- In lossless mode, raising `--credit` or lowering `--inbox` moves more of the Router's pool into paths that can stall; a multi-node deadlock was found and fixed that way once (`FINDINGS.md`, session 8). Retest multi-node after touching `--n-recv`, `--credit` or `--inbox`.

## Workflow

Solo research repo (single remote, single author) — no PR process. Commit directly to feature branches (e.g. `omp_reboot`); keep it informal.
