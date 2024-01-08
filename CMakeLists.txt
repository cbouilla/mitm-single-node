include(CMakePrintHelpers) # Print variables easily

cmake_minimum_required(VERSION 3.20)


project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json
# sorry, for now gcc is the source of truth.


# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)
add_compile_options(-march=native -ggdb -Wall -Wextra -pedantic)
# Needed this, other wise it won't link
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")


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


add_executable(speck32_demo
               demos/speck32_demo.cpp
	       mitm_sequential.hpp
	       include/dict.hpp
               demos/speck32.hpp)

add_executable(sha2_demo
        demos/sha2_demo.cpp
	demos/sha256.c
        mitm_sequential.hpp
        include/dict.hpp)

add_executable(sha2-x86_demo
        demos/sha2_demo.cpp
	demos/sha256_x86.c
        mitm_sequential.hpp
        include/dict.hpp)


# # Link TBB to your executable
# target_link_libraries(speck32_demo PRIVATE TBB::tbb)
# #[[
# add_executable(test demos/basic_demo.cpp mitm_sequential.hpp
#         include/dict.hpp
#         include/dict.hpp
#         demos/basic_demo.cpp)]]
