###### SHA256

add_executable(sha2_collision_demo
    sha2_collision_demo.cpp sha256.c)
target_include_directories(sha2_collision_demo PRIVATE ../include)

add_executable(sha2_claw_demo
    sha2_claw_demo.cpp sha256.c)
target_include_directories(sha2_claw_demo PRIVATE ../include)


###### Speck

add_executable(double_speck64_demo 
    double_speck64_demo.cpp)
target_include_directories(double_speck64_demo PRIVATE ../include)

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



###### DES

add_executable(double_DES_demo usuba_des.c double_DES_demo.cpp)
target_include_directories(double_DES_demo PRIVATE ../include)
target_link_libraries(double_DES_demo PRIVATE OpenSSL::Crypto)


add_executable(mpi_double_DES_bench usuba_des.c mpi_double_DES_bench.cpp)
target_include_directories(mpi_double_DES_bench PRIVATE ../include)
target_link_libraries(mpi_double_DES_bench PRIVATE OpenSSL::Crypto)
target_link_libraries(mpi_double_DES_bench PUBLIC MPI::MPI_CXX)


###### AES

add_executable(double_aes_demo 
    aes.c double_aes_demo.cpp)
target_include_directories(double_aes_demo PRIVATE ../include)

add_executable(mpi_double_aes_bench
    aes.c mpi_double_aes_bench.cpp)
target_include_directories(mpi_double_aes_bench PRIVATE ../include)
target_link_libraries(mpi_double_aes_bench PUBLIC MPI::MPI_CXX)

