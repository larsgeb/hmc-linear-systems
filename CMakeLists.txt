cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)
#set (CMAKE_CXX_FLAGS "-O4")

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SOURCE_FILES generateSynthetics.cpp linearalgebra.cpp auxiliary.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp)
add_executable(generateSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp)
add_executable(runSampling ${SOURCE_FILES})

set(SOURCE_FILES generateGradient.cpp linearalgebra.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp)
add_executable(generateGradient ${SOURCE_FILES})

set(SOURCE_FILES testMultivariate.cpp linearalgebra.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/SparseVector.cpp
        SparseLinearAlgebra/src/AlgebraLib/SparseMatrix.cpp SparseLinearAlgebra/src/AlgebraLib/Matrix.cpp
        SparseLinearAlgebra/src/AlgebraLib/Matrix.hpp SparseLinearAlgebra/src/AlgebraLib/Vector.cpp
        SparseLinearAlgebra/src/AlgebraLib/Vector.hpp SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp
        SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.cpp SparseLinearAlgebra/src/AlgebraLib/FullAlgebra.hpp
        SparseLinearAlgebra/src/AlgebraLib/Globals.hpp)
add_executable(testMultivariate ${SOURCE_FILES})
