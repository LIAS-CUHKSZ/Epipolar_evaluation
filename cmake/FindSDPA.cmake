# https://github.com/coin-or/Gravity
# https://github.com/coin-or/Gravity/blob/master/cmake/FindSDPA.cmake

# ATTENTION: CHANGE THIS TO YOUR OWN PATH
set(SDPA_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/sdpa-7.3.18")
message("Looking for Sdpa in ${SDPA_ROOT_DIR}")

find_path(SDPA_INCLUDE_DIR	sdpa_call.h	HINTS "${SDPA_ROOT_DIR}")
find_path(Mumps_INCLUDE_DIR dmumps_c.h HINTS "${SDPA_ROOT_DIR}/mumps/build/include")
find_library (SDPA_LIBRARY 	libsdpa.a HINTS "${SDPA_ROOT_DIR}/")
find_library(Mumps_LIBRARY1 libdmumps.a	 HINTS "${SDPA_ROOT_DIR}/mumps/build/lib" )
find_library(Mumps_LIBRARY2 libmumps_common.a  HINTS "${SDPA_ROOT_DIR}/mumps/build/lib" )
find_library(Mumps_LIBRARY3 libpord.a  HINTS "${SDPA_ROOT_DIR}/mumps/build/lib" )
find_library(Mumps_LIBRARY4 libmpiseq.a	 HINTS "${SDPA_ROOT_DIR}/mumps/build/libseq" )
find_library(BLAS_LIBRARY libopenblas.a HINTS  "/usr/lib/x86_64-linux-gnu/")
find_library(FORTRAN_LIBRARY libgfortran.so.5 HINTS "/usr/lib/x86_64-linux-gnu/")
find_library(FORTRAN_LIBRARY2 libquadmath.so.0	HINTS "/usr/lib/x86_64-linux-gnu/")


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SDPA DEFAULT_MSG SDPA_LIBRARY SDPA_INCLUDE_DIR)

if(SDPA_FOUND)
	message("—- Found Sdpa under ${SDPA_INCLUDE_DIR}")
    set(SDPA_INCLUDE_DIRS ${SDPA_INCLUDE_DIR}  ${Mumps_INCLUDE_DIR})
    set(SDPA_LIBRARIES ${SDPA_LIBRARY} ${Mumps_LIBRARY1} ${Mumps_LIBRARY2}
        ${Mumps_LIBRARY3} ${Mumps_LIBRARY4} ${BLAS_LIBRARY} "pthread" ${FORTRAN_LIBRARY} ${FORTRAN_LIBRARY2})
    message("—- Set Sdpa lib  ${SDPA_LIBRARY}")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(SDPA_LIBRARIES "${SDPA_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(SDPA_FOUND)

mark_as_advanced(SDPA_LIBRARY SDPA_INCLUDE_DIR)
