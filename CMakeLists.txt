cmake_minimum_required(VERSION 3.20)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json

# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)

add_compile_options(-march=native -ggdb -Wall)
# The little bird says clang has a nicer error report than gcc!

# Only the MPI C API is used, from C++ (no MPI-2 C++ bindings, which MPI-3 deprecated and
# Open MPI 5 removed).  Asking for the C component alone spares us FindMPI's mpicxx probe,
# which fails on clusters whose C++ wrapper is missing or broken.
find_package(MPI REQUIRED COMPONENTS C)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)     # DES implementation

# hwloc 2.x: NUMA placement of the threads (include/tools.hpp).  Hard requirement, no
# fallback.  On a cluster: module load hwloc, or HWLOC_ROOT=/path cmake ...
find_path(HWLOC_INCLUDE_DIR hwloc.h REQUIRED HINTS ENV HWLOC_ROOT PATH_SUFFIXES include)
find_library(HWLOC_LIBRARY hwloc REQUIRED HINTS ENV HWLOC_ROOT PATH_SUFFIXES lib lib64)
add_library(hwloc::hwloc UNKNOWN IMPORTED)
set_target_properties(hwloc::hwloc PROPERTIES IMPORTED_LOCATION "${HWLOC_LIBRARY}"
                      INTERFACE_INCLUDE_DIRECTORIES "${HWLOC_INCLUDE_DIR}")
message(STATUS "hwloc: ${HWLOC_LIBRARY}")

# The router tests run under ctest (examples/CMakeLists.txt).  Open MPI binds a rank to one core at
# np <= 2, which serialises a spinning OpenMP team: unbind by default there.
enable_testing()
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version OUTPUT_VARIABLE MPIEXEC_VERSION_OUTPUT
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if(MPIEXEC_VERSION_OUTPUT MATCHES "Open MPI" AND NOT MPIEXEC_PREFLAGS)
    set(MPIEXEC_PREFLAGS "--bind-to;none" CACHE STRING "mpiexec flags before the program" FORCE)
endif()

add_subdirectory(examples)

include("config/avx2.cmake")
include("config/avx512.cmake")