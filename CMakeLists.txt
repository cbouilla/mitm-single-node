cmake_minimum_required(VERSION 3.20)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json

# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)

add_compile_options(-march=native -ggdb -Wall)
# The little bird says clang has a nicer error report than gcc!

find_package(OpenMP REQUIRED)

#find_package(MPI REQUIRED)
#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} #ARGS")
#
## Check if MPI was found
#if(MPI_CXX_FOUND)
#  # Link the MPI libraries to your targets
#  target_link_libraries(sha2_collision_demo PUBLIC MPI::MPI_CXX)
#  target_link_libraries(sha2_claw_demo PUBLIC MPI::MPI_CXX)
#endif()
#
## link OpenMP (node level parallelism) and MPI (multiple nodes parallelism)
#target_link_libraries(sha2_collision_demo PUBLIC MPI::MPI_CXX)
#
#target_link_libraries(sha2_collision_demo PUBLIC MPI::MPI_CXX)

add_subdirectory(examples)