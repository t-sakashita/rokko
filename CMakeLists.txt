#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2018 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

#
# Basic setup
#

cmake_minimum_required(VERSION 3.1)
project(rokko NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

if(CONFIG)
  message(STATUS "Loading configration: " ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

#
# Options
#

include(components)
components_init()

# RPATH setting
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

# C++, C, and Fortran compilers
enable_language(CXX C)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
if(BUILD_Fortran_BINDING)
  enable_language(Fortran)
  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
  endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # ancient gfortran does not support iso_c_bindings
  try_compile(CMake_Fortran_SUPPORTS_ISO_C_BINDING
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/config/check_iso_c_binding.f90
    OUTPUT_VARIABLE LOG)
  if(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
    set(BUILD_Fortran_BINDING OFF)
    MESSAGE(STATUS "Do not build Fortran binding since the fortran compiler does not support iso_c_binding.")
  endif(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
endif(BUILD_Fortran_BINDING)

# Linker for Fortran Program
if(BUILD_Fortran_BINDING)
  message(STATUS "Use C++ linker for Fortran program: " ${USE_CXX_LINKER})
  set(ROKKO_USE_CXX_LINKER ${USE_CXX_LINKER})
endif(BUILD_Fortran_BINDING)

# OpenMP
set(ROKKO_HAVE_OPENMP False)
if(NOT DISABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(ROKKO_HAVE_OPENMP True)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Almost always OpenMP flags are same both for C and for Fortran.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif(OPENMP_FOUND)
endif(NOT DISABLE_OPENMP)

# Python
if(BUILD_Python_BINDING)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})
endif(BUILD_Python_BINDING)

# lapack library
if((EXISTS $ENV{MKLROOT}) AND (NOT BLA_VENDOR))
  set(BLA_VENDOR Intel10_64lp)
endif((EXISTS $ENV{MKLROOT}) AND (NOT BLA_VENDOR))
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK not found")
endif(NOT LAPACK_FOUND)

# cblas library
try_compile(_HAVE_CBLAS
  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/config/check_cblas.c
  LINK_LIBRARIES ${BLAS_LIBRARIES} OUTPUT_VARIABLE LOG)
if(_HAVE_CBLAS)
  message(STATUS "Found CBLAS in BLAS")
  set(BUILD_CBLAS OFF)
else(_HAVE_CBLAS)
  try_compile(_HAVE_CBLAS_MKL ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/config/check_cblas_mkl.c
    LINK_LIBRARIES ${BLAS_LIBRARIES} OUTPUT_VARIABLE LOG)
  if(_HAVE_CBLAS_MKL)
    message(STATUS "Found CBLAS in MKL")
    set(BUILD_CBLAS OFF)
    set(HAVE_CBLAS_MKL TRUE)
  endif(_HAVE_CBLAS_MKL)
endif(_HAVE_CBLAS)
if(BUILD_CBLAS)
  set(CBLAS_LIBRARIES cblas)
else(BUILD_CBLAS)
  set(CBLAS_LIBRARIES "")
endif(BUILD_CBLAS)

# lapacke library
find_library(_LAPACKE_LIBRARY NAMES lapacke DOC "The LAPACKE library")
if(_LAPACKE_LIBRARY)
  message(STATUS "LAPACKE library: ${_LAPACKE_LIBRARY}")
  list(APPEND LAPACK_LIBRARIES ${_LAPACKE_LIBRARY})
endif(_LAPACKE_LIBRARY)
try_compile(_HAVE_LAPACKE ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/config/check_lapacke.c
  LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} OUTPUT_VARIABLE LOG)
if(_HAVE_LAPACKE)
  message(STATUS "Found LAPACKE")
  set(BUILD_LAPACKE OFF)
else(_HAVE_LAPACKE)
  try_compile(_HAVE_LAPACKE_MKL ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/config/check_lapacke_mkl.c
    LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} OUTPUT_VARIABLE LOG)
  if(_HAVE_LAPACKE_MKL)
    message(STATUS "Found LAPACKE in MKL")
    set(BUILD_LAPACKE OFF)
    set(HAVE_LAPACKE_MKL TRUE)
  endif(_HAVE_LAPACKE_MKL)
endif(_HAVE_LAPACKE)
if(BUILD_LAPACKE)
  set(LAPACKE_LIBRARIES lapacke)
else(BUILD_LAPACKE)
  set(LAPACKE_LIBRARIES "")
endif(BUILD_LAPACKE)

# boost Library
set(ROKKO_HAVE_BOOST False)
if(BUILD_Python_BINDING)
  find_package(Boost COMPONENTS python system)
else(BUILD_Python_BINDING)
  find_package(Boost COMPONENTS system)
endif(BUILD_Python_BINDING)
if(Boost_INCLUDE_DIR)
  set(ROKKO_HAVE_BOOST True)
  include_directories(BEFORE ${Boost_INCLUDE_DIR})
endif(Boost_INCLUDE_DIR)
components_check()

# MPI library
if(ROKKO_HAVE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
  else(MPI_FOUND)
    set(ROKKO_HAVE_MPI False)
  endif(MPI_FOUND)
endif(ROKKO_HAVE_MPI)
components_check()

#
# parallel dense solvers
#

# find ScaLAPACK
if(ROKKO_HAVE_SCALAPACK)
  find_package(ScaLAPACK)
  set(ROKKO_HAVE_SCALAPACK ${SCALAPACK_FOUND})
  string(REGEX REPLACE ";" " " SCAPALACK_LIBS "${SCAPALACK_LIBRARIES}")
endif(ROKKO_HAVE_SCALAPACK)
components_check()

# find EigenExa
if(ROKKO_HAVE_EIGENEXA)
  find_package(EigenExa)
  set(ROKKO_HAVE_EIGENEXA ${EIGENEXA_FOUND})
  include_directories(${EIGENEXA_INCLUDE_DIR})
  string(REGEX REPLACE ";" " " EIGENEXA_LIBS "${EIGENEXA_LIBRARIES}")
endif(ROKKO_HAVE_EIGENEXA)
components_check()

# find Elemental
if(ROKKO_HAVE_ELEMENTAL)
  find_package(Elemental)
  set(ROKKO_HAVE_ELEMENTAL ${ELEMENTAL_FOUND})
  include_directories(${ELEMENTAL_INCLUDE_DIR})
  string(REGEX REPLACE ";" " " ELEMENTAL_LIBS "${ELEMENTAL_LIBRARIES}")
endif(ROKKO_HAVE_ELEMENTAL)
if(ROKKO_HAVE__ELEMENTAL)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  message(STATUS "Compiler option: ${CMAKE_CXX_FLAGS}")
endif(ROKKO_HAVE__ELEMENTAL)
components_check()

# find ELPA
if(ROKKO_HAVE_ELPA)
  find_package(ELPA)
  set(ROKKO_HAVE_ELPA ${ELPA_FOUND})
  string(REGEX REPLACE ";" " " ELPA_LIBS "${ELPA_LIBRARIES}")
endif(ROKKO_HAVE_ELPA)
components_check()

#
# parallel sparse solvers
#

# find PETSc/SLEPc
if(ROKKO_HAVE_PETSC)
  find_package(PETSc COMPONENTS CXX)
  set(ROKKO_HAVE_PETSC ${PETSC_FOUND})
  include_directories(${PETSC_INCLUDE_DIR})
  string(REGEX REPLACE ";" " " PETSC_LIBS "${PETSC_LIBRARIES}")
endif(ROKKO_HAVE_PETSC)
if(ROKKO_HAVE_SLEPC)
  find_package(SLEPc COMPONENTS CXX)
  set(ROKKO_HAVE_SLEPC ${SLEPC_FOUND})
  include_directories(${SLEPC_INCLUDE_DIR})
  string(REGEX REPLACE ";" " " SLEPC_LIBS "${SLEPC_LIBRARIES}")
endif(ROKKO_HAVE_SLEPC)
components_check()

# find Anasazi
if(ROKKO_HAVE_ANASAZI)
  find_package(Anasazi)
  set(ROKKO_HAVE_ANASAZI ${ANASAZI_FOUND})
  include_directories(${ANASAZI_INCLUDE_DIR})
  string(REGEX REPLACE ";" " " ANASAZI_LIBS "${ANASAZI_LIBRARIES}")
endif(ROKKO_HAVE_ANASAZI)
components_check()

set(ROKKO_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")

configure_file(${PROJECT_SOURCE_DIR}/rokko/config.h.in
  ${PROJECT_BINARY_DIR}/rokko/config.h)
if(NOT EXISTS ${PROJECT_BINARY_DIR}/rokko/mangling.h)
  include(FortranCInterface)
  FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/rokko/mangling.h
                           MACRO_NAMESPACE "ROKKO_"
                           SYMBOL_NAMESPACE "ROKKO_")
  if(NOT FortranCInterface_COMPILED)
    configure_file(${PROJECT_SOURCE_DIR}/rokko/mangling.h.in
    ${PROJECT_BINARY_DIR}/rokko/mangling.h COPYONLY)
  endif(NOT FortranCInterface_COMPILED)
endif(NOT EXISTS ${PROJECT_BINARY_DIR}/rokko/mangling.h)
install(FILES ${PROJECT_BINARY_DIR}/rokko/config.h ${PROJECT_BINARY_DIR}/rokko/mangling.h
  COMPONENT headers DESTINATION include/rokko/)
install(FILES cmatrix.h DESTINATION include COMPONENT headers)

configure_file(${PROJECT_SOURCE_DIR}/config/rokko-config.cmake.in
  ${PROJECT_BINARY_DIR}/config/rokko-config.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/config/include.mk.in
  ${PROJECT_BINARY_DIR}/config/include.mk @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/config/rokko-config.cmake
  ${PROJECT_BINARY_DIR}/config/include.mk
  ${PROJECT_SOURCE_DIR}/config/use-rokko.cmake
  ${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake
  ${PROJECT_SOURCE_DIR}/config/run_test.cmake
  ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake DESTINATION share/rokko/)

# standard include directories
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
if(BUILD_CBLAS)
  include_directories(${PROJECT_BINARY_DIR}/3rd-party/cblas/include)
  include_directories(${PROJECT_SOURCE_DIR}/3rd-party/cblas/include)
else(BUILD_CBLAS)
  if(HAVE_CBLAS_MKL)
    include_directories(${PROJECT_SOURCE_DIR}/3rd-party/cblas_mkl)
  endif(HAVE_CBLAS_MKL)
endif(BUILD_CBLAS)
if(BUILD_LAPACKE)
  include_directories(${PROJECT_BINARY_DIR}/3rd-party/lapacke/include)
  include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke/include)
else(BUILD_LAPACKE)
  if(HAVE_LAPACKE_MKL)
    include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke_mkl)
  endif(HAVE_LAPACKE_MKL)
endif(BUILD_LAPACKE)
if(BUILD_BASE_CXX)
  include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)
endif(BUILD_BASE_CXX)

# test environment
if(BUILD_TESTS OR BUILD_EXAMPLES)
  set(GOOGLETEST_VERSION 1)
  add_subdirectory(3rd-party/googletest)
  include_directories(3rd-party/googletest/include)
  include(add_rokko_test)
  enable_testing()
  include(CTest)
endif(BUILD_TESTS OR BUILD_EXAMPLES)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(tool)
if(BUILD_TESTS)
  add_subdirectory(test)
endif(BUILD_TESTS)
if(BUILD_EXAMPLES)
  add_subdirectory(example)
endif(BUILD_EXAMPLES)
if(BUILD_TUTORIALS)
  add_subdirectory(tutorial)
endif(BUILD_TUTORIALS)
if(BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif(BUILD_BENCHMARKS)
if(BUILD_DEVELOPERS)
  add_subdirectory(developer)
endif(BUILD_DEVELOPERS)

components_summary()
