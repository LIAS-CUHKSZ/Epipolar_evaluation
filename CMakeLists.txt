project(EPP)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3  REQUIRED)
find_package(OpenCV  REQUIRED)
find_package(SDPA  REQUIRED)

add_definitions(-DM_PI=3.1415926535f)

include_directories( ${SDPA_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} )
include_directories("D:/Desktop/epipolar_eval/sdpa/mumps/build/include") # path to your mumps include folder
# by default, when building SDPA, the mumps lib will be built in the sdpa/mumps/build/lib folder

# Settings.
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build type not specified, using RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_definitions("-O3 -msse2 -msse3 -march=native -std=c++17")
set( CMAKE_CXX_FLAGS "-march=native -fopenmp" )
include_directories(
    ./src
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -funroll-all-loops")
set(CMAKE_Fortran_FLAGS "-Wall -fPIC -funroll-all-loops")

# Executable.
add_executable(epipolar_eval
    src/cameras.cc
    src/images.cc
    src/eval.cc
    src/pt3d.hpp
    src/const_est.hpp
    src/manifold_GN.hpp
    src/utils.hpp
    src/npt_pose.cpp
    src/construct_coeff_constraint.cpp
    src/construct_coeff_objective.cpp
)

target_link_libraries(epipolar_eval ${OpenCV_LIBS} ${SDPA_LIBRARIES} )

# m-consistent evaluation 
add_executable(MonteCarlo src/m_pts_eval.cc
src/cameras.cc
    src/images.cc
    src/pt3d.hpp
    src/const_est.hpp
    src/manifold_GN.hpp
    src/utils.hpp
    src/npt_pose.cpp
    src/construct_coeff_constraint.cpp
    src/construct_coeff_objective.cpp)
  target_link_libraries(MonteCarlo ${OpenCV_LIBS} ${SDPA_LIBRARIES} )
