# Sequential meet in the middle

```
                AbstractDomain.hpp
                        +
                       / \
                      /   \
                     /     \
                    /       \
                   /         \
AbstractClawProblem.hpp     AbstractCollisionProblem.hpp
                  ^           ^
                  |           |
         claw_egine.hpp    collision_engine.hpp
                  ^           ^
                   \         /
                    engine.hpp
                        ^
                        |
                     mitm.hpp
                        ^
                        |
                   your_code.cpp
                   /    |    \
follows AbstractDomain and {AbstractClawProblem or AbstractCollisionProblem})
```
