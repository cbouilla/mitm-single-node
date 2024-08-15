
add_executable(sha2_collision_demo
    sha2_collision_demo.cpp
	sha256_x86.c)
target_link_libraries(sha2_collision_demo PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(sha2_collision_demo PRIVATE ../include)

add_executable(sha2_claw_demo
    sha2_claw_demo.cpp
	sha256_x86.c)
target_link_libraries(sha2_claw_demo PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(sha2_claw_demo PRIVATE ../include)

add_executable(double_speck32_demo
    double_speck32_demo.cpp
    speck3264.c)
target_link_libraries(double_speck32_demo PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(double_speck32_demo PRIVATE ../include)

add_executable(double_speck64_demo
    double_speck64_demo.cpp
    speck64128.c)
target_link_libraries(double_speck64_demo PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(double_speck64_demo PRIVATE ../include)

if(MPI_CXX_FOUND)
    add_executable(mpi_double_speck64_demo
        mpi_double_speck64_demo.cpp
        speck64128.c)
    target_link_libraries(mpi_double_speck64_demo PUBLIC OpenMP::OpenMP_CXX)
    target_include_directories(mpi_double_speck64_demo PRIVATE ../include)
    target_link_libraries(mpi_double_speck64_demo PUBLIC MPI::MPI_CXX)

    add_executable(mpi_naive_double_speck64_demo
        mpi_naive_double_speck64_demo.cpp
        speck64128.c)
    target_include_directories(mpi_naive_double_speck64_demo PRIVATE ../include)
    target_link_libraries(mpi_naive_double_speck64_demo PUBLIC MPI::MPI_CXX)
endif()
