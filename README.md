**Note** Developement for parallel `mitm` moved to [mitm](https://github.com/akaalharbi/mitm). This repository will be freezed to use it as a correctness test against the parallel version.

# Sequential meet in the middle

```
                AbstractDomain.hpp
                        |
                        x
                       / \
                      /   \
                     /     \
                    /       \
                   /         \
                  /           \
                 /             \
	        /               \
AbstractClawProblem.hpp     AbstractCollisionProblem.hpp
              /                  /
             /                  /
            /                  / 
           / naive_engien.hpp /
          /     |            /
         /      |           /
claw_egine.hpp  |    collision_engine.hpp 
         \      |         /
          \     |        /
           \    |       /
            \   |      /
prng.hpp  <- engine.hpp -> dict.hpp
                |
              mitm.hpp
                |
            your_code.cpp
 

```
