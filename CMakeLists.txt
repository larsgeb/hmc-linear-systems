cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)
set (CMAKE_CXX_FLAGS "-O4")

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SOURCE_FILES generateSynthetics.cpp
        AlgebraLib/src/algebra_lib/sparse_algebra.cpp AlgebraLib/src/algebra_lib/sparse_vector.cpp
        AlgebraLib/src/algebra_lib/sparse_matrix.cpp AlgebraLib/src/algebra_lib/matrix.cpp
        AlgebraLib/src/algebra_lib/matrix.hpp AlgebraLib/src/algebra_lib/vector.cpp
        AlgebraLib/src/algebra_lib/vector.hpp AlgebraLib/src/algebra_lib/algebra_lib.hpp
        AlgebraLib/src/algebra_lib/full_algebra.cpp AlgebraLib/src/algebra_lib/full_algebra.hpp
        AlgebraLib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)
add_executable(generateSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp src/random/randomnumbers.cpp
        AlgebraLib/src/algebra_lib/sparse_algebra.cpp AlgebraLib/src/algebra_lib/sparse_vector.cpp
        AlgebraLib/src/algebra_lib/sparse_matrix.cpp AlgebraLib/src/algebra_lib/matrix.cpp
        AlgebraLib/src/algebra_lib/matrix.hpp AlgebraLib/src/algebra_lib/vector.cpp
        AlgebraLib/src/algebra_lib/vector.hpp AlgebraLib/src/algebra_lib/algebra_lib.hpp
        AlgebraLib/src/algebra_lib/full_algebra.cpp AlgebraLib/src/algebra_lib/full_algebra.hpp
        AlgebraLib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)
add_executable(runSampling ${SOURCE_FILES})

set(SOURCE_FILES generateGradient.cpp src/random/randomnumbers.cpp
        AlgebraLib/src/algebra_lib/sparse_algebra.cpp AlgebraLib/src/algebra_lib/sparse_vector.cpp
        AlgebraLib/src/algebra_lib/sparse_matrix.cpp AlgebraLib/src/algebra_lib/matrix.cpp
        AlgebraLib/src/algebra_lib/matrix.hpp AlgebraLib/src/algebra_lib/vector.cpp
        AlgebraLib/src/algebra_lib/vector.hpp AlgebraLib/src/algebra_lib/algebra_lib.hpp
        AlgebraLib/src/algebra_lib/full_algebra.cpp AlgebraLib/src/algebra_lib/full_algebra.hpp
        AlgebraLib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)
add_executable(generateGradient ${SOURCE_FILES})

set(SOURCE_FILES testMultivariate.cpp src/random/randomnumbers.cpp
        AlgebraLib/src/algebra_lib/sparse_algebra.cpp AlgebraLib/src/algebra_lib/sparse_vector.cpp
        AlgebraLib/src/algebra_lib/sparse_matrix.cpp AlgebraLib/src/algebra_lib/matrix.cpp
        AlgebraLib/src/algebra_lib/matrix.hpp AlgebraLib/src/algebra_lib/vector.cpp
        AlgebraLib/src/algebra_lib/vector.hpp AlgebraLib/src/algebra_lib/algebra_lib.hpp
        AlgebraLib/src/algebra_lib/full_algebra.cpp AlgebraLib/src/algebra_lib/full_algebra.hpp
        AlgebraLib/src/algebra_lib/globals.hpp src/hmc/hmc.hpp src/hmc/prior.cpp src/hmc/prior.hpp src/hmc/data.cpp src/hmc/data.hpp src/hmc/posterior.cpp src/hmc/posterior.hpp src/hmc/sampler.cpp src/hmc/sampler.hpp src/hmc/forward_model.cpp src/hmc/forward_model.hpp)
add_executable(testMultivariate ${SOURCE_FILES})
