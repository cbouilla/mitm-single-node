include(CMakePrintHelpers) # Print variables easily

cmake_minimum_required(VERSION 3.20)


project(mitm_sequential)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json

# for AES header
set(AES_IMPL "aesni" CACHE STRING "Choose an AES Implementation")
set_property(CACHE AES_IMPL PROPERTY STRINGS auto portable aesni neon)


# sorry, for now gcc is the source of truth.


# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)
# from -g -no-pie etc, are flags for boost::backtrace
add_compile_options(-march=native -ggdb -Wall -Wextra -pedantic -O0  -no-pie -fno-pie -rdynamic -DBOOST_STACKTRACE_USE_ADDR2LINE)
# Needed this, other wise it won't link
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
# Using almost every possible sanitizer!
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb -fsanitize=address,undefined,leak,pointer-compare,pointer-subtract,shift,integer-divide-by-zero,null,return,vla-bound,signed-integer-overflow,bounds,alignment,bool,enum,vptr")


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




add_executable(sha2_collision_demo
        examples/sha2_collision_demo.cpp
	examples/sha256_x86.c # examples/sha256.c
	mitm.hpp)

target_compile_definitions(sha2_collision_demo PUBLIC COLLISION_DEBUG)


add_executable(sha2_claw_demo
        examples/sha2_claw_demo.cpp
	examples/sha256_x86.c
	mitm.hpp)

# Debugging the execution if the results don't add up
target_compile_definitions(sha2_claw_demo PUBLIC CLAW_DEBUG)



find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  # target_link_libraries(sha2_collision_1collision_demo PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(sha2_collision_demo PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(sha2_claw_demo PUBLIC OpenMP::OpenMP_CXX)
endif()



find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
# Check if MPI was found
if(MPI_CXX_FOUND)
  # Link the MPI libraries to your targets
  target_link_libraries(sha2_collision_demo PUBLIC MPI::MPI_CXX)
  target_link_libraries(sha2_claw_demo PUBLIC MPI::MPI_CXX)
endif()

# link OpenMP (node level parallelism) and MPI (multiple nodes parallelism)
target_link_libraries(sha2_collision_demo
                      PUBLIC MPI::MPI_CXX)

target_link_libraries(sha2_collision_demo
                      PUBLIC MPI::MPI_CXX)
# # Link TBB to your executable

# #[[
# add_executable(test examples/basic_demo.cpp mitm_sequential.hpp
#         dict.hpp
#         dict.hpp
#         examples/basic_demo.cpp)]]


# add this block to after `add_executable(...)`
if("${AES_IMPL}" STREQUAL "aesni")
    target_compile_definitions(sha2_claw_demo PUBLIC USE_INTEL_AESNI)
    if(MSVC)
        target_compile_options(sha2_claw_demo PRIVATE /arch:SSE2)
    else()
        target_compile_options(sha2_claw_demo PRIVATE -maes)
    endif()
elseif("${AES_IMPL}" STREQUAL "neon")
    target_compile_definitions(sha2_claw_demo PUBLIC USE_NEON_AES)
    target_compile_options(sha2_claw_demo PRIVATE -march=armv8-a+crypto)
elseif("${AES_IMPL}" STREQUAL "portable")
    target_compile_definitions(sha2_claw_demo PUBLIC USE_CXX_AES)
endif()


if("${AES_IMPL}" STREQUAL "aesni")
    target_compile_definitions(sha2_collision_demo PUBLIC USE_INTEL_AESNI)
    if(MSVC)
        target_compile_options(sha2_collision_demo PRIVATE /arch:SSE2)
    else()
        target_compile_options(sha2_collision_demo PRIVATE -maes)
    endif()
elseif("${AES_IMPL}" STREQUAL "neon")
    target_compile_definitions(sha2_collision_demo PUBLIC USE_NEON_AES)
    target_compile_options(sha2_collision_demo PRIVATE -march=armv8-a+crypto)
elseif("${AES_IMPL}" STREQUAL "portable")
    target_compile_definitions(sha2_collision_demo PUBLIC USE_CXX_AES)
endif()
 
