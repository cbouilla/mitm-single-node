cmake_minimum_required(VERSION 3.20)
project(mitm_sequential)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Compiler command json
# sorry, for now gcc is the source of truth.
add_compile_options(-ggdb -Wall -Wextra -pedantic)


include_directories(include)


add_executable(demo_speck32
        demos/speck32_demo.cpp
        mitm_sequential.hpp
        include/dict.hpp
        demos/speck32.hpp)

#[[
add_executable(test demos/basic_demo.cpp mitm_sequential.hpp
        include/dict.hpp
        include/dict.hpp
        demos/basic_demo.cpp)]]
