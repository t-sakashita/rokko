#
# This module is provided as ROKKO_USE_FILE by RokkoConfig.cmake.  It can
# be included in a project to load the needed compiler and linker
# settings to use Rokko.
#

if(NOT ROKKO_USE_FILE_INCLUDED)
  set(ROKKO_USE_FILE_INCLUDED 1)
  list(APPEND CMAKE_MODULE_PATH ${ROKKO_ROOT_DIR}/share/rokko)

  # compilers and common options
  if(NOT PREVENT_ROKKO_COMPILERS)
    set(CMAKE_C_COMPILER ${ROKKO_CMAKE_C_COMPILER} CACHE FILEPATH "C compiler." FORCE)
    set(CMAKE_C_FLAGS ${ROKKO_CMAKE_C_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS_DEBUG ${ROKKO_CMAKE_C_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_C_FLAGS_RELEASE ${ROKKO_CMAKE_C_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_CXX_COMPILER ${ROKKO_CMAKE_CXX_COMPILER} CACHE FILEPATH "CXX compiler." FORCE)
    set(CMAKE_CXX_FLAGS ${ROKKO_CMAKE_CXX_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_CXX_FLAGS_DEBUG ${ROKKO_CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE ${ROKKO_CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_Fortran_COMPILER ${ROKKO_CMAKE_Fortran_COMPILER} CACHE FILEPATH "Fortran compiler." FORCE)
    set(CMAKE_Fortran_FLAGS ${ROKKO_CMAKE_Fortran_FLAGS} CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG ${ROKKO_CMAKE_Fortran_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE ${ROKKO_CMAKE_Fortran_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during release builds." FORCE)
  endif(NOT PREVENT_ROKKO_COMPILERS)

  # OpenMP
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  endif(OPENMP_FOUND)

  # MPI library
  find_package(MPI)
  set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
  set(ROKKO_LIBRARIES ${ROKKO_LIBRARIES} ${MPI_CXX_LIBRARIES})

  # Add include directories needed to use ROKKO and dependent libraries
  include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${ROKKO_INCLUDE_DIR} ${ROKKO_Boost_INCLUDE_DIR} ${ROKKO_EIGENEXA_INCLUDE_DIR} ${ROKKO_ELEMENTAL_INCLUDE_DIR} ${ROKKO_ANASAZI_INCLUDE_DIR} ${ROKKO_PETSC_INCLUDE_DIR} ${ROKKO_SLEPC_INCLUDE_DIR})

  # Add link directories needed to use ROKKO and dependent libraries
  link_directories(${ROKKO_ROOT_DIR}/lib)

  # RPATH setting
  set(CMAKE_INSTALL_NAME_DIR "${ROKKO_ROOT_DIR}/lib" FORCE)
  set(CMAKE_INSTALL_RPATH "${ROKKO_ROOT_DIR}/lib" FORCE)

  # test macro
  include(${ROKKO_ROOT_DIR}/share/rokko/add_rokko_test.cmake)
endif(NOT ROKKO_USE_FILE_INCLUDED)
