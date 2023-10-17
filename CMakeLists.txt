cmake_minimum_required(VERSION 3.26)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)

# add_executable(mitm_sequential demos/main.cpp mitm_sequential.hpp)
add_executable(mitm_sequential demos/demo.cpp mitm_sequential.hpp)
add_executable(crax_mitm demos/crax-demo.cpp mitm_sequential.hpp)

