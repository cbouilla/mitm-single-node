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


# NOT BUILT: naive_double_speck64_demo.cpp, on top of include/naive/.  That is the
# naive all-to-all MITM, a different engine; it has not been ported to the one-rank-
# per-node topology and does not compile.  Kept as the starting point for that work.
