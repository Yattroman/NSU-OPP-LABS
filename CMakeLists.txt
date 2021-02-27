cmake_minimum_required(VERSION 3.7)
project(OPP_NSU_LABS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")

include_directories(/usr/include/mpi/)

set(LAB1V1_SOURCE laboratory-1/version-1/main.cpp laboratory-1/version-1/mv/mvInit.cpp laboratory-1/version-1/mv/mvOperations.cpp)
set(LAB1V2_SOURCE laboratory-1/version-2/main.cpp)
set(LAB1V3_SOURCE laboratory-1/version-2/main.cpp)

add_executable(laboratory-1-v1 ${LAB1V1_SOURCE} )
add_executable(laboratory-1-v2 ${LAB1V2_SOURCE} )
add_executable(laboratory-1-v3 ${LAB1V2_SOURCE} )