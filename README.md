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
                  /           \
AbstractClawProblem.hpp     AbstractCollisionProblem.hpp
                  \           /
                   \         /
                    \       /
                     \     /
                      \   /
                       \ /
                        ^
                        |
      prng.hpp  <- base_engine.hpp -> dict.hpp
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
                   your_code.cpp
 
follows AbstractDomain and {AbstractClawProblem or AbstractCollisionProblem})
```
