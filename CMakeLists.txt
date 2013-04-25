cmake_minimum_required(VERSION 2.8)

project(rokko NONE)
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

# find Elemental
#set(ELEMENTAL_DIR "/home/sakashita/elemental-0.78_installed")
#set(ELEMENTAL_DIR "/opt/nano/rokko")
find_package(Elemental)
message(STATUS  "Display ELEM_INC:  ${ELEM_INC}")
message(STATUS  "Display ELEM_LIB:  ${ELEM_LIB}")
message(STATUS  "Display ELEM_LIBS:  ${ELEM_LIBS}")

# directories
add_subdirectory(sample_scalapack)
add_subdirectory(test)
add_subdirectory(sample_eigenK)
add_subdirectory(sample_elemental)




