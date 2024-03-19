**Note** Developement for parallel `mitm` moved to [mitm](https://github.com/akaalharbi/mitm). 

**Note**: 
To compile:
```
rm -rf CMakeCache.txt  CMakeFiles/ && cmake -S . -B . -D AES_IMPL=aesni && make && ./sha2_claw_demo
```

**Note**: This repository is kept seperated from the parallel implementation to continue testing on various `mitm` scenarios where we need a working version. Most development will be in the `demos` folder.

# Sequential meet in the middle

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
