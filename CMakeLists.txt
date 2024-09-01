cmake_minimum_required(VERSION 3.20)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json

# common compiler option (between clang and gcc, not sure about the others)
# add_compile_options(-ggdb -Wall -Wextra -pedantic -fsanitize=address)

add_compile_options(-march=native -ggdb -Wall)
# The little bird says clang has a nicer error report than gcc!

find_package(MPI REQUIRED)
find_package(OpenSSL REQUIRED)     # DES implementation


add_subdirectory(examples)

include("config/avx2.cmake")
include("config/avx512.cmake")