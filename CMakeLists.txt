cmake_minimum_required(VERSION 3.16)
project(EPP)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3  REQUIRED)
find_package(OpenCV  REQUIRED)
find_package(SDPA  REQUIRED)
set(opengv_DIR ${CMAKE_SOURCE_DIR}/src/3rdparty/opengv/build)
find_package(opengv REQUIRED)
set(OpenGV_LIBS ${CMAKE_SOURCE_DIR}/src/3rdparty/opengv/build/lib/libopengv.so)
find_package(yaml-cpp REQUIRED)
find_package(OpenMP REQUIRED)


# add_definitions(-DREMOVE_OUTLIER)

include_directories( ${SDPA_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIR} /home/neo/Epipolar_evaluation/src/3rdparty/opengv/include )


if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build type not specified, using RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_definitions("-O3 -msse2 -msse3 -march=native -std=c++17")
set( CMAKE_CXX_FLAGS "-march=native -fopenmp" )
include_directories(
    ./src
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -funroll-all-loops ${OpenMP_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "-Wall -fPIC -funroll-all-loops ${OpenMP_Fortran_FLAGS}")


include_directories(src/data_loader src/methods src/methods/E_sdp )

file(GLOB_RECURSE SRCLIST 
  src/data_loader/*.cpp 
  src/data_loader/*.cc
  src/data_loader/*.hpp
  src/methods/*.cpp 
  src/methods/*.cc 
  src/methods/*.hpp
  src/utils/hpp
)

# Create a library from SRCLIST
add_library(eval_lib STATIC ${SRCLIST})
target_link_libraries(eval_lib ${OpenCV_LIBS} ${SDPA_LIBRARIES} ${OpenGV_LIBS} ${YAML_CPP_LIBRARIES} OpenMP::OpenMP_CXX)

# Executable, eval for epipolar geometry estimation
# add_executable(epipolar_eval src/eval.cc src/utils.hpp)
# target_link_libraries(epipolar_eval eval_lib)


# add_executable(kitti_eval src/kitti_eval.cc src/utils.hpp)
# target_link_libraries(kitti_eval eval_lib)

# add_executable(eth3dGN src/eth3d_GN.cc src/utils.hpp)
# target_link_libraries(eth3dGN eval_lib)

# add_executable(kitti_analysis src/kitti_analysis.cc src/utils.hpp)
# target_link_libraries(kitti_analysis eval_lib)

add_executable(MonteCarlo src/m_pts_eval.cc src/utils.hpp)
target_link_libraries(MonteCarlo eval_lib)
