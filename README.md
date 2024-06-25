

**Note**: 
To compile:
```
rm -rf CMakeCache.txt  CMakeFiles/ && cmake -S . -B . -D AES_IMPL=aesni && make && ./sha2_claw_demo
```

**Note**: Parallel Meet-in-The-Middle (mitm). For a single process mitm, please visit the branch `single-node`


#  meet in the middle

```
                AbstractDomain.hpp
                        ^
                        |
                        x
                       / \
                      /   \
                     /     \
                    /       \
                   /         \
AbstractClawProblem.hpp     AbstractCollisionProblem.hpp
                   \         /
                    \       /
                     \     /
                      \   /
                       \ /
                        ^
                        |
        prng.hpp  <- engine.hpp -> dict.hpp
                        ^
                        |
                        x
                       / \
                      /   \
                     |     |
                  +---      --+
                  |           |
         claw_egine.hpp    collision_engine.hpp
                  ^           ^
                   \         /
                    --+   +--
                       \ /
    		        |       common_mpi.hpp
                        x        /
                       / \      /
                      /   \    /
                     |     |  /
                 +----     - / -+
                 |          /   |
                 |         /    |
                 |        /     |
          receiver.hpp-> + <-sender.hpp
                  ^           ^
                   \         /
                    --+   +--
                       \ /
                        |
                parallel_engine.hpp
                        |
                     mitm.hpp
                        ^
                        |
                   your_code.cpp

```



```bash
├── mitm.hpp
├── docs (Follow this for your implementation)
│   ├── AbstractClawProblem.hpp
│   ├── AbstractCollisionProblem.hpp
│   └── AbstractDomain.hpp
├── include
│   ├── claw_engine.hpp
│   ├── collision_engine.hpp
│   ├── common.hpp
│   ├── counters.hpp
│   ├── dict.hpp
│   ├── engine.hpp
│   ├── hash_table
│   ├── mpi_common.hpp
│   ├── naive_engine.hpp
│   ├── parallel_engine.hpp
│   ├── permutations (you can use examples here for mixing function)
│   │   └── AES.hpp
│   ├── receiver.hpp
│   ├── sender.hpp
│   └── util (independent code used in mitm)
│       ├── folder_creator.hpp
│       ├── memory.hpp
│       ├── prng.hpp
│       └── timing.hpp
├── CMakeLists.txt
├── examples (that uses claw and collision finding)
│   ├── sha256.c
│   ├── sha256_x86.c
│   ├── sha2_claw_demo.cpp
│   ├── sha2_claw_demo.cpp
│   ├── sha2_collision_demo.cpp
│   ├── bits_lib.hpp
│   ├── scripts (ignore this)
│   └── └── ...
├── playground (ignore this)
│   ├── cachelines.cpp
│   ├── collisions_summary.csv
│   ├── collision_summary.csv
│   ├── process_collisions_summary.py
│   ├── summary_collision.csv
│   └── summary_collisions.ipynb
├── data (ignore this)
│   └── ...
└── README.md
```
