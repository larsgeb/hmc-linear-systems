cmake_minimum_required(VERSION 2.8)
project(hmc-linear-system)
#set (CMAKE_CXX_FLAGS "-O4")

set (CMAKE_CXX_FLAGS "-std=gnu++11 -g ${CMAKE_CXX_FLAGS}")

SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SOURCE_FILES generateSynthetics.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp src/HMC/HMC.hpp src/HMC/Prior.cpp src/HMC/Prior.hpp src/HMC/Data.cpp src/HMC/Data.hpp src/HMC/Posterior.cpp src/HMC/Posterior.hpp src/HMC/Sampler.cpp src/HMC/Sampler.hpp src/HMC/ForwardModel.cpp src/HMC/ForwardModel.hpp)
add_executable(generateSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp randomnumbers.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp src/HMC/HMC.hpp src/HMC/Prior.cpp src/HMC/Prior.hpp src/HMC/Data.cpp src/HMC/Data.hpp src/HMC/Posterior.cpp src/HMC/Posterior.hpp src/HMC/Sampler.cpp src/HMC/Sampler.hpp src/HMC/ForwardModel.cpp src/HMC/ForwardModel.hpp)
add_executable(runSampling ${SOURCE_FILES})

set(SOURCE_FILES generateGradient.cpp randomnumbers.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp src/HMC/HMC.hpp src/HMC/Prior.cpp src/HMC/Prior.hpp src/HMC/Data.cpp src/HMC/Data.hpp src/HMC/Posterior.cpp src/HMC/Posterior.hpp src/HMC/Sampler.cpp src/HMC/Sampler.hpp src/HMC/ForwardModel.cpp src/HMC/ForwardModel.hpp)
add_executable(generateGradient ${SOURCE_FILES})

set(SOURCE_FILES testMultivariate.cpp randomnumbers.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp src/HMC/HMC.hpp src/HMC/Prior.cpp src/HMC/Prior.hpp src/HMC/Data.cpp src/HMC/Data.hpp src/HMC/Posterior.cpp src/HMC/Posterior.hpp src/HMC/Sampler.cpp src/HMC/Sampler.hpp src/HMC/ForwardModel.cpp src/HMC/ForwardModel.hpp)
add_executable(testMultivariate ${SOURCE_FILES})
