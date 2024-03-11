include(CMakePrintHelpers) # Print variables easily

cmake_minimum_required(VERSION 3.20)


project(mitm_sequential)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json
# sorry, for now gcc is the source of truth.


# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)
# from -g -no-pie etc, are flags for boost::backtrace
add_compile_options(-march=native -ggdb -Wall -Wextra -pedantic -O0 -g -no-pie -fno-pie -rdynamic -DBOOST_STACKTRACE_USE_ADDR2LINE)
# Needed this, other wise it won't link
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb -fsanitize=address")


# Boost stacktrace # sudo apt-get install libboost-stacktrace-dev
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -no-pie -fno-pie -lboost_stacktrace_backtrace -DBOOST_STACKTRACE_LINK -ldl -DBOOST_STACKTRACE_USE_ADDR2LINE -D BOOST_STACKTRACE_LINK")

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


# add_executable(sha2_collision_1collision_demo
#         demos/sha2_collision_1collision_demo.cpp
# 	demos/sha256_x86.c # demos/sha256.c
# 	mitm.hpp
# 	engine.hpp
# 	naive_engine.hpp
# 	collision_engine.hpp
# 	AbstractCollisionProblem.hpp
# 	AbstractDomain.hpp
#         dict.hpp
#       	include/memory.hpp
# 	include/prng.hpp
# 	include/timing.hpp)



# add_executable(sha2_collision_demo
#         demos/sha2_collision_demo.cpp
# 	demos/sha256_x86.c # demos/sha256.c
# 	mitm.hpp
# 	engine.hpp
# 	naive_engine.hpp
# 	collision_engine.hpp
# 	AbstractCollisionProblem.hpp
# 	AbstractDomain.hpp
#         dict.hpp
#       	include/memory.hpp
# 	include/prng.hpp
# 	include/timing.hpp)



add_executable(sha2_claw_demo
        demos/sha2_claw_demo.cpp
	demos/sha256_x86.c
	mitm.hpp
	engine.hpp
	naive_engine.hpp
	collision_engine.hpp
	AbstractCollisionProblem.hpp
	AbstractDomain.hpp
        dict.hpp
      	include/memory.hpp
	include/prng.hpp
	include/timing.hpp)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  # target_link_libraries(sha2_collision_1collision_demo PUBLIC OpenMP::OpenMP_CXX)
  # target_link_libraries(sha2_collision_demo PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(sha2_claw_demo PUBLIC OpenMP::OpenMP_CXX)
endif()


# # Link TBB to your executable

# #[[
# add_executable(test demos/basic_demo.cpp mitm_sequential.hpp
#         dict.hpp
#         dict.hpp
#         demos/basic_demo.cpp)]]
