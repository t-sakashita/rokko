cmake_minimum_required(VERSION 2.8)

project(frank NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# enable C and C++ compilers
#enable_language(CXX)
enable_language(C CXX Fortran)

find_package(MPI)

# find Eigen3 Library
find_package(Eigen3)

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(BEFORE ${EIGEN3_INCLUDE_DIR})

set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})

# directories
add_subdirectory(sample_scalapack)
add_subdirectory(test)
add_subdirectory(sample_eigenK)




