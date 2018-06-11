cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lopenblas -fopenmp")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic") # Want to debug?
#set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Werror") # Want to debug really bad?

include_directories(../armadillo-code/include) # or whatever your current Armadillo directory is

set(SOURCE_FILES_SAMPLER src/executables/runSampling.cpp src/random/randomnumbers.cpp src/hmc/linearSampler.cpp src/hmc/linearSampler.hpp)
set(SOURCE_FILES_QUADRATIC src/executables/createQuadraticForm.cpp src/random/randomnumbers.cpp src/hmc/linearSampler.cpp src/hmc/linearSampler.hpp)

add_executable(hmc_sampler ${SOURCE_FILES_SAMPLER})
add_executable(quadratic ${SOURCE_FILES_QUADRATIC})

target_link_libraries(hmc_sampler openblas)
target_link_libraries(quadratic openblas)
