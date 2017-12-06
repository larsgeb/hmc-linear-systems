cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)
#set (CMAKE_CXX_FLAGS "-O2")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE release)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_C_FLAGS "fopenmp")

find_package(LAPACK)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES src/hmc/runSampling.cpp src/random/randomnumbers.cpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp
        src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp
        src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)
add_executable(runSampling ${SOURCE_FILES})

target_link_libraries(runSampling ${lapackblas_libraries} ${ARMADILLO_LIBARIES})
target_link_libraries(runSampling ${ARMADILLO_LIBRARIES})
