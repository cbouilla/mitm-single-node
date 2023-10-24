cmake_minimum_required(VERSION 3.14)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)
include_directories(include)

add_executable(test demos/basic_demo.cpp mitm_sequential.hpp
        include/dict.hpp
        include/dict.hpp
        demos/basic_demo.cpp)