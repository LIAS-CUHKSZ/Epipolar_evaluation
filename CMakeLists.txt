project(ETH3DFormatLoader)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

# External dependency: Eigen3
find_package(Eigen3  REQUIRED)
find_package(OpenCV  REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} )

# Settings.
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build type not specified, using RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_definitions("-O3 -msse2 -msse3 -std=c++17")
include_directories(
    ./src
)

# Executable.
add_executable(ETH3DFormatLoader
    src/cameras.cc
    src/images.cc
    src/main.cc
    src/pt3d.hpp
    src/const_est.hpp
    src/manifold_GN.hpp
    src/utils.hpp
)

target_link_libraries(ETH3DFormatLoader ${OpenCV_LIBS} )
