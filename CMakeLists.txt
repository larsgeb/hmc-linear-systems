cmake_minimum_required(VERSION 3.7)
project(hmc_vsp)

set(CMAKE_CXX_STANDARD 11)

#set(SOURCE_FILES sampling.cpp aux.hpp aux.cpp montecarlo.hpp montecarlo.cpp auxiliary.cpp auxiliary.hpp linearalgebra.cpp linearalgebra.hpp randomnumbers.cpp randomnumbers.hpp samplingNew.cpp)
#add_executable(sampling ${SOURCE_FILES})

set(SOURCE_FILES makeSynthetics.cpp linearalgebra.cpp toyForwardModel.cpp auxiliary.cpp)
add_executable(makeSynthetics ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp linearalgebra.cpp randomnumbers.cpp toyForwardModel.cpp auxiliary.cpp montecarlo.cpp)
add_executable(runSamplingToyModel ${SOURCE_FILES})

set(SOURCE_FILES runSampling.cpp linearalgebra.cpp randomnumbers.cpp tomographyForwardModel.cpp auxiliary.cpp montecarlo.cpp)
add_executable(runSamplingTomography ${SOURCE_FILES})
