cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lopenblas -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic") # Want to debug?
#set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Werror") # Want to debug really bad?

include_directories(../armadillo-code/include) # Or whatever your current Armadillo directory is

set(SOURCE_FILES src/hmc/runSampling.cpp src/random/randomnumbers.cpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp
        src/hmc/data.cpp src/hmc/data.hpp src/hmc/sampler.cpp
        src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)

add_executable(hmc_sampler ${SOURCE_FILES})
target_link_libraries(hmc_sampler openblas)
