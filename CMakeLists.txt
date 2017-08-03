cmake_minimum_required(VERSION 3.7)
project(hmc-linear-system)
#set (CMAKE_CXX_FLAGS "-O4")

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES makeSynthetics.cpp linearalgebra.cpp auxiliary.cpp)
add_executable(makeSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp linearalgebra.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp)
add_executable(runSampling ${SOURCE_FILES})

set(SOURCE_FILES generateGradient.cpp linearalgebra.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp)
add_executable(generateGradient ${SOURCE_FILES})

set(SOURCE_FILES testMultivariate.cpp linearalgebra.cpp randomnumbers.cpp auxiliary.cpp montecarlo.cpp)
add_executable(testMultivariate ${SOURCE_FILES})
