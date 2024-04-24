

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
                        |
                     mitm.hpp
                        ^
                        |
                   your_code.cpp

```
