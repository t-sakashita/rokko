cmake_minimum_required(VERSION 2.8)

project(rokko NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(add_rokko_test)

# enable C and C++ compilers
enable_language(C CXX Fortran)


# options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_EIGEN_S "Build eigen_s library" ON)
option(BUILD_EIGEN_SX "Build eigen_sx library" ON)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # Almost always OpenMP flags are same both for C and for Fortran.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# MPI library
find_package(MPI)

# Eigen3 Library
find_package(Eigen3)

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(BEFORE ${EIGEN3_INCLUDE_DIR})
if (BOOST_INCLUDE_DIR)
  include_directories(BEFORE ${BOOST_INCLUDE_DIR})
endif()

set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})

# find Elemental
#set(ELEMENTAL_DIR "/home/sakashita/elemental-0.78_installed")
#set(ELEMENTAL_DIR "/opt/nano/rokko")
find_package(Elemental)
message(STATUS  "Display ELEM_INC:  ${ELEM_INC}")
message(STATUS  "Display ELEM_LIB:  ${ELEM_LIB}")
message(STATUS  "Display ELEM_LIBS:  ${ELEM_LIBS}")

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(sample_scalapack)
# add_subdirectory(test)
add_subdirectory(sample_eigenK)
add_subdirectory(sample_elemental)
add_subdirectory(benchmark)
