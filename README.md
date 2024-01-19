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
 
follows AbstractDomain and {AbstractClawProblem or AbstractCollisionProblem})
```
