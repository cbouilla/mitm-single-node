include(CMakePrintHelpers) # Print variables easily

cmake_minimum_required(VERSION 3.20)


project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json
# sorry, for now gcc is the source of truth.


# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)
add_compile_options(-march=native -ggdb -Wall -Wextra -pedantic -O0 -ggdb )
# Needed this, other wise it won't link
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb -fsanitize=address")


# The little bird says clang has a nicer error report than gcc!
# export CC=/usr/bin/clang
# export CXX=/usr/bin/clang++
# set(CMAKE_C_COMPILER /usr/bin/clang)
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)


cmake_print_variables(CMAKE_CXX_COMPILER_ID)

# it's a bit ugly, since it won't generate agnostic Makefile
# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#     # GCC specific options
#     add_compile_options( -fanalyzer)
# elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     # Clang specific options
#     add_compile_options(-Rpass=.*)
# endif()

# Find the TBB package
# find_package(TBB REQUIRED)

include_directories(include)


# add_executable(speck32_demo
#                demos/speck32_demo.cpp
# 	       claw_engine.hpp
#        	       base_engine.hpp
# 	       AbstractClawProblem.hpp
# 	       AbstractDomain.hpp
# 	       dict.hpp
# 	       base_engine.hpp
# 	       include/memory.hpp
# 	       include/prng.hpp
# 	       include/timing.hpp
#                demos/speck32.hpp)

add_executable(sha2_collision_demo
        demos/sha2_collision_demo.cpp
	demos/sha256_x86.c # demos/sha256.c
	naive_engine.hpp
	collision_engine.hpp
	base_engine.hpp
	AbstractCollisionProblem.hpp
	AbstractDomain.hpp
        dict.hpp
      	include/memory.hpp
	include/prng.hpp
	include/timing.hpp)



add_executable(sha2_claw_demo
        demos/sha2_claw_demo.cpp
	demos/sha256_x86.c
	naive_engine.hpp
	collision_engine.hpp
	base_engine.hpp
	AbstractCollisionProblem.hpp
	AbstractDomain.hpp
        dict.hpp
      	include/memory.hpp
	include/prng.hpp
	include/timing.hpp)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(sha2_collision_demo PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(sha2_claw_demo PUBLIC OpenMP::OpenMP_CXX)
endif()


# # Link TBB to your executable

# #[[
# add_executable(test demos/basic_demo.cpp mitm_sequential.hpp
#         dict.hpp
#         dict.hpp
#         demos/basic_demo.cpp)]]
