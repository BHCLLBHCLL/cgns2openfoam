cmake_minimum_required(VERSION 3.14)
project(cgns2openfoam VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default paths from environment or user-specified
set(HDF5_ROOT "$ENV{HDF5_ROOT}" CACHE PATH "HDF5 installation root")
set(CGNS_ROOT "$ENV{CGNS_ROOT}" CACHE PATH "CGNS installation root")

# Fallback paths if environment not set
if(NOT HDF5_ROOT)
    set(HDF5_ROOT "/home/sdcll221/Programs/hdf5-2.0.0")
endif()
if(NOT CGNS_ROOT)
    set(CGNS_ROOT "/home/sdcll221/Programs/CGNS-4.5.1")
endif()

message(STATUS "HDF5_ROOT: ${HDF5_ROOT}")
message(STATUS "CGNS_ROOT: ${CGNS_ROOT}")

# HDF5
set(HDF5_INCLUDE_DIR "${HDF5_ROOT}/include")
set(HDF5_LIBRARY "${HDF5_ROOT}/lib/libhdf5.so")

if(NOT EXISTS "${HDF5_INCLUDE_DIR}/hdf5.h")
    message(FATAL_ERROR "HDF5 header not found at ${HDF5_INCLUDE_DIR}")
endif()
if(NOT EXISTS "${HDF5_LIBRARY}")
    message(FATAL_ERROR "HDF5 library not found at ${HDF5_LIBRARY}")
endif()

# CGNS
set(CGNS_INCLUDE_DIR "${CGNS_ROOT}/include")
set(CGNS_LIBRARY "${CGNS_ROOT}/lib/libcgns.so")

if(NOT EXISTS "${CGNS_INCLUDE_DIR}/cgnslib.h")
    message(FATAL_ERROR "CGNS header not found at ${CGNS_INCLUDE_DIR}")
endif()
if(NOT EXISTS "${CGNS_LIBRARY}")
    message(FATAL_ERROR "CGNS library not found at ${CGNS_LIBRARY}")
endif()

message(STATUS "HDF5 include: ${HDF5_INCLUDE_DIR}")
message(STATUS "HDF5 library: ${HDF5_LIBRARY}")
message(STATUS "CGNS include: ${CGNS_INCLUDE_DIR}")
message(STATUS "CGNS library: ${CGNS_LIBRARY}")

# Executable
add_executable(cgns2openfoam
    src/main.cpp
    src/cgns_reader.cpp
    src/openfoam_writer.cpp
)

target_include_directories(cgns2openfoam PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CGNS_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
)

target_link_libraries(cgns2openfoam PRIVATE
    ${CGNS_LIBRARY}
    ${HDF5_LIBRARY}
)

# RPATH for finding shared libraries at runtime
set_target_properties(cgns2openfoam PROPERTIES
    INSTALL_RPATH "${CGNS_ROOT}/lib;${HDF5_ROOT}/lib"
    BUILD_RPATH "${CGNS_ROOT}/lib;${HDF5_ROOT}/lib"
)

# Install
install(TARGETS cgns2openfoam DESTINATION bin)
