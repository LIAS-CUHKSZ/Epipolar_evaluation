
cmake_minimum_required(VERSION 3.16)
project(EPP)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3  REQUIRED)
find_package(OpenCV  REQUIRED)
find_package(SDPA  REQUIRED)

# you can also add config.cmake of opengv to camke prefix and use find_package(OpenGV REQUIRED)
set(OpenGV_LIBS "D:/Desktop/epipolar_eval/opengv/lib/libopengv.a")
include_directories("D:/Desktop/epipolar_eval/opengv/include")

# path to your mumps include folder, for windows user only
# by default, when building SDPA, the mumps lib will be built in the sdpa/mumps/build/lib folder
# if you install SDPA from software repo, just comment this line
include_directories("D:/Desktop/epipolar_eval/sdpa/mumps/build/include") 

include_directories( ${SDPA_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} )


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


include_directories(src/data_loader src/methods src/methods/E_sdp)

file(GLOB_RECURSE SRCLIST 
  src/data_loader/*.cpp 
  src/data_loader/*.cc
  src/data_loader/*.hpp
  src/methods/*.cpp 
  src/methods/*.cc 
  src/methods/*.hpp
  src/utils/hpp
)

# Executable,eval for epipolar geometry estimation
add_executable(epipolar_eval src/eval.cc
    src/utils.hpp
    ${SRCLIST}
)
target_link_libraries(epipolar_eval ${OpenCV_LIBS} ${SDPA_LIBRARIES} ${OpenGV_LIBS} )

# m-consistent evaluation 
# add_executable(MonteCarlo src/m_pts_eval.cc
#                src/utils.hpp
#                ${SRCLIST}
# )
# target_link_libraries(MonteCarlo ${OpenCV_LIBS} ${SDPA_LIBRARIES}  ${OpenGV_LIBS})
