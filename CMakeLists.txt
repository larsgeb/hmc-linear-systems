cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)
#set (CMAKE_CXX_FLAGS "-O4")

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SOURCE_FILES generateSynthetics.cpp
        algebra_lib/src/algebra_lib/sparse_algebra.cpp algebra_lib/src/algebra_lib/sparse_vector.cpp
        algebra_lib/src/algebra_lib/sparse_matrix.cpp algebra_lib/src/algebra_lib/matrix.cpp
        algebra_lib/src/algebra_lib/matrix.hpp algebra_lib/src/algebra_lib/vector.cpp
        algebra_lib/src/algebra_lib/vector.hpp algebra_lib/src/algebra_lib/algebra_lib.hpp
        algebra_lib/src/algebra_lib/full_algebra.cpp algebra_lib/src/algebra_lib/full_algebra.hpp
        algebra_lib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forwardmodel.cpp src/hmc/forwardmodel.hpp)
add_executable(generateSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp randomnumbers.cpp
        algebra_lib/src/algebra_lib/sparse_algebra.cpp algebra_lib/src/algebra_lib/sparse_vector.cpp
        algebra_lib/src/algebra_lib/sparse_matrix.cpp algebra_lib/src/algebra_lib/matrix.cpp
        algebra_lib/src/algebra_lib/matrix.hpp algebra_lib/src/algebra_lib/vector.cpp
        algebra_lib/src/algebra_lib/vector.hpp algebra_lib/src/algebra_lib/algebra_lib.hpp
        algebra_lib/src/algebra_lib/full_algebra.cpp algebra_lib/src/algebra_lib/full_algebra.hpp
        algebra_lib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forwardmodel.cpp src/hmc/forwardmodel.hpp)
add_executable(runSampling ${SOURCE_FILES})

set(SOURCE_FILES generateGradient.cpp randomnumbers.cpp
        algebra_lib/src/algebra_lib/sparse_algebra.cpp algebra_lib/src/algebra_lib/sparse_vector.cpp
        algebra_lib/src/algebra_lib/sparse_matrix.cpp algebra_lib/src/algebra_lib/matrix.cpp
        algebra_lib/src/algebra_lib/matrix.hpp algebra_lib/src/algebra_lib/vector.cpp
        algebra_lib/src/algebra_lib/vector.hpp algebra_lib/src/algebra_lib/algebra_lib.hpp
        algebra_lib/src/algebra_lib/full_algebra.cpp algebra_lib/src/algebra_lib/full_algebra.hpp
        algebra_lib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forwardmodel.cpp src/hmc/forwardmodel.hpp)
add_executable(generateGradient ${SOURCE_FILES})

set(SOURCE_FILES testMultivariate.cpp randomnumbers.cpp
        algebra_lib/src/algebra_lib/sparse_algebra.cpp algebra_lib/src/algebra_lib/sparse_vector.cpp
        algebra_lib/src/algebra_lib/sparse_matrix.cpp algebra_lib/src/algebra_lib/matrix.cpp
        algebra_lib/src/algebra_lib/matrix.hpp algebra_lib/src/algebra_lib/vector.cpp
        algebra_lib/src/algebra_lib/vector.hpp algebra_lib/src/algebra_lib/algebra_lib.hpp
        algebra_lib/src/algebra_lib/full_algebra.cpp algebra_lib/src/algebra_lib/full_algebra.hpp
        algebra_lib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forwardmodel.cpp src/hmc/forwardmodel.hpp)
add_executable(testMultivariate ${SOURCE_FILES})
