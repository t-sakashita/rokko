cmake_minimum_required(VERSION 2.8)

project(frank NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# enable C and C++ compilers
enable_language(CXX)

find_package(MPI)

# find Eigen3 Library
find_package(Eigen3)

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(BEFORE ${EIGEN3_INCLUDE_DIR})

# directories
add_subdirectory(sample_scalapack)




