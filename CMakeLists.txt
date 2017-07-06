cmake_minimum_required(VERSION 3.7)
project(hmc_vsp)

set(CMAKE_CXX_STANDARD 11)

set(EXE1 make_data.cpp aux.hpp aux.cpp)
add_executable(make_data ${EXE1})

set(SOURCE_FILES sampling.cpp aux.hpp aux.cpp mc.hpp mc.cpp)
add_executable(sampling ${SOURCE_FILES})