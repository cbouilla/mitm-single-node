# Every driver runs the same engine (MPI + OpenMP), so they all link the same way.
function(mitm_example name)
    add_executable(${name} ${name}.cpp ${ARGN})
    target_include_directories(${name} PRIVATE ../include)
    target_link_libraries(${name} PRIVATE MPI::MPI_C OpenMP::OpenMP_CXX Threads::Threads hwloc::hwloc)
    # Keep <mpi.h> from pulling in the deprecated MPI-2 C++ bindings; we never use them.
    target_compile_definitions(${name} PRIVATE OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX _MPICC_H)
endfunction()


###### Speck

mitm_example(double_speck64_demo)
mitm_example(double_speck64_bench)


###### DES

mitm_example(double_DES_demo usuba_des.cpp)
target_link_libraries(double_DES_demo PRIVATE OpenSSL::Crypto)

mitm_example(double_DES_bench usuba_des.cpp)
target_link_libraries(double_DES_bench PRIVATE OpenSSL::Crypto)


###### AES

mitm_example(double_aes_demo aes.c)
mitm_example(double_aes_bench aes.c)


###### SHA256

mitm_example(sha2_claw_demo sha256.c)
mitm_example(sha2_collision_demo sha256.c)



###### The Router: its test suite and its benchmark

mitm_example(router_test)
target_compile_definitions(router_test PRIVATE ROUTER_PARANOID)
mitm_example(router_bench)

# one rank, two ranks, four ranks: 8 spinning threads in all on an 8-core laptop (12 at np 4, the
# minimum of one service thread, one sender and one receiver per rank).  np 1 pins (one rank owns the
# machine); np 2 and np 4 put several ranks on the one host, whose masks overlap, so they pass --no-bind
# (the Router refuses to pin overlapping masks) -- pinned multi-rank needs --map-by numa on a real cluster.
add_test(NAME router_np1 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
         $<TARGET_FILE:router_test> --senders 4 --receivers 3 --quiet)
add_test(NAME router_np2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
         $<TARGET_FILE:router_test> --senders 2 --receivers 1 --no-bind --quiet)
add_test(NAME router_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         $<TARGET_FILE:router_test> --senders 1 --receivers 1 --no-bind --quiet)
