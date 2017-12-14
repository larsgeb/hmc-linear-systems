cmake_minimum_required(VERSION 3.8)
project(hmc-linear-system)
set (CMAKE_CXX_FLAGS "-O2")

set(CMAKE_CXX_FLAGS "-std=gnu++11 -g ${CMAKE_CXX_FLAGS}")

SET(CMAKE_BUILD_TYPE debug)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_C_FLAGS "fopenmp -lopenblas -llapack -std=c++0x")

#find_package(Armadillo REQUIRED)
include_directories(/cluster/home/lgebraad/armadillo/armadillo-8.200.2/include/armadillo_bits)
include_directories(/cluster/home/lgebraad/armadillo/armadillo-8.200.2/include)
include_directories(/cluster/home/lgebraad/armadillo/armadillo-8.200.2/)

set(SOURCE_FILES src/hmc/runSampling.cpp src/random/randomnumbers.cpp src/hmc/hmc.hpp src/hmc/prior.cpp
        src/hmc/prior.hpp
        src/hmc/data.cpp src/hmc/data.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forward_model.cpp
        src/hmc/forward_model.hpp)
add_executable(hmc_sampler ${SOURCE_FILES})

find_package(OpenBLAS)
#find_package(Armadillo 2.4.2 REQUIRED)

target_link_libraries(hmc_sampler ${CMAKE_CURRENT_SOURCE_DIR}/armadillo/armadillo-8.200.2/libarmadillo.so)
target_link_libraries(hmc_sampler openblas)

