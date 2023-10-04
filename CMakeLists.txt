cmake_minimum_required(VERSION 3.14)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)

add_executable(test demo.cpp mitm_sequential.hpp)