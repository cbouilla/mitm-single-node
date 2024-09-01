add_executable(sha2_collision_demo
    sha2_collision_demo.cpp sha256.c)
target_include_directories(sha2_collision_demo PRIVATE ../include)

add_executable(sha2_claw_demo
    sha2_claw_demo.cpp sha256.c)
target_include_directories(sha2_claw_demo PRIVATE ../include)

add_executable(double_speck64_demo 
    double_speck64_demo.cpp)
target_include_directories(double_speck64_demo PRIVATE ../include)


add_executable(double_des_demo double_des_demo.cpp)
target_include_directories(double_des_demo PRIVATE ../include des)
target_link_libraries(double_des_demo PRIVATE OpenSSL::Crypto)


add_executable(mpi_double_speck64_demo
    mpi_double_speck64_demo.cpp)
target_include_directories(mpi_double_speck64_demo PRIVATE ../include)
target_link_libraries(mpi_double_speck64_demo PUBLIC MPI::MPI_CXX)

add_executable(mpi_naive_double_speck64_demo
    mpi_naive_double_speck64_demo.cpp)
target_include_directories(mpi_naive_double_speck64_demo PRIVATE ../include)
target_link_libraries(mpi_naive_double_speck64_demo PUBLIC MPI::MPI_CXX)

add_executable(mpi_double_speck64_bench
    mpi_double_speck64_bench.cpp)
target_include_directories(mpi_double_speck64_bench PRIVATE ../include)
target_link_libraries(mpi_double_speck64_bench PUBLIC MPI::MPI_CXX)

add_executable(mpi_sha2_big_collision_demo
    mpi_sha2_big_collision_demo.cpp sha256.c)
target_include_directories(mpi_sha2_big_collision_demo PRIVATE ../include)
target_link_libraries(mpi_sha2_big_collision_demo PUBLIC MPI::MPI_CXX)
