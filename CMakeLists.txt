cmake_minimum_required(VERSION 3.7)
project(OPP_NSU_LABS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")

include_directories(/usr/include/mpi/)

set(LAB1V1_SOURCE laboratory-1/version-1/main.cpp laboratory-1/version-1/mv/mvInit_l1_1st.cpp laboratory-1/version-1/mv/mvOperations_l1_1st.cpp)
set(LAB1V2_SOURCE laboratory-1/version-2/main.cpp laboratory-1/version-2/mv/mvInit_l1_2nd.cpp laboratory-1/version-2/mv/mvOperations_l1_2nd.cpp)
set(LAB1V3_SOURCE laboratory-1/version-3/main.cpp laboratory-1/version-3/mv/mvInit_l1_3rd.cpp laboratory-1/version-3/mv/mvOperations_l1_3rd.cpp)

set(LAB2V1_SOURCE laboratory-2/version-1/main.cpp laboratory-2/version-1/mv/mvInit_l2_1st.cpp laboratory-2/version-1/mv/mvOperations_l2_1st.cpp)
set(LAB2V2_SOURCE laboratory-2/version-2/main.cpp laboratory-2/version-2/mv/mvInit_l2_2nd.cpp laboratory-2/version-2/mv/mvOperations_l2_2nd.cpp)
set(LAB2V3_SOURCE laboratory-2/version-3/main.cpp laboratory-2/version-3/mv/mvInit_l2_3rd.cpp laboratory-2/version-3/mv/mvOperations_l2_3rd.cpp)

set(LAB3_SOURCE laboratory-3/main_l3.cpp)

set(LAB4_SOURCE laboratory-4/main_l4.cpp)

set(LAB5_SOURCE laboratory-5/lab5.cpp laboratory-5/main.cpp)

add_executable(laboratory-3 ${LAB3_SOURCE})
add_executable(laboratory-4 ${LAB4_SOURCE})
add_executable(laboratory-5 ${LAB5_SOURCE})