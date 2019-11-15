#
# This module is provided as ROKKO_USE_FILE by rokko-config.cmake.  It can
# be included in a project to load the needed compiler and linker
# settings to use Rokko.
#

if(NOT ROKKO_USE_FILE_INCLUDED)
  set(ROKKO_USE_FILE_INCLUDED 1)
  list(APPEND CMAKE_MODULE_PATH ${ROKKO_ROOT_DIR}/share/rokko)

  # compilers and common options
  if(NOT PREVENT_ROKKO_COMPILERS)
    set(CMAKE_BUILD_TYPE ${ROKKO_CMAKE_BUILD_TYPE} CACHE STRING "Type of build" FORCE)
    message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

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
    set(CMAKE_BUILD_TYPE ${ROKKO_CMAKE_BUILD_TYPE} CACHE STRING "Type of build" FORCE)

    if(ROKKO_Fortran_BINDING)
      enable_language(CXX C Fortran)
    else(ROKKO_Fortran_BINDING)
      enable_language(CXX C)
    endif(ROKKO_Fortran_BINDING)

    # OpenMP
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif(OPENMP_FOUND)

    # MPI library
    find_package(MPI)
    set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
    set(ROKKO_LIBRARIES ${ROKKO_LIBRARIES} ${MPI_CXX_LIBRARIES})

  endif(NOT PREVENT_ROKKO_COMPILERS)

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

macro(display_rokko_components)
  message(STATUS "Boost_INCLUDE_DIR=${ROKKO_Boost_INCLUDE_DIR}")

  message(STATUS "Solvers:")
  # serial dense solvers
  message(STATUS "  BLAS_LIBRARIES=${ROKKO_BLAS_LIBRARIES}")
  message(STATUS "  LAPACK_LIBRARIES=${ROKKO_LAPACK_LIBRARIES}")

  # parallel dense solvers
  message(STATUS "  ROKKO_HAVE_SCALAPACK:" ${ROKKO_HAVE_SCALAPACK})
  if(ROKKO_HAVE_SCALAPACK)
    message(STATUS "  SCALAPACK_LIBRARIES=${ROKKO_SCALAPACK_LIBRARIES}")
  endif(ROKKO_HAVE_SCALAPACK)

  message(STATUS "  ROKKO_HAVE_EIGENEXA:" ${ROKKO_HAVE_EIGENEXA})
  if(ROKKO_HAVE_EIGENEXA)
    message(STATUS "  EIGENEXA_INCLUDE_DIR=${ROKKO_EIGENEXA_INCLUDE_DIR}")
    message(STATUS "  EIGENEXA_LIBRARIES=${ROKKO_EIGENEXA_LIBRARIES}")
  endif(ROKKO_HAVE_EIGENEXA)

  message(STATUS "  ROKKO_HAVE_ELPA:" ${ROKKO_HAVE_ELPA})
  if(ROKKO_HAVE_ELPA)
    message(STATUS "  ELPA_INCLUDE_DIR=${ROKKO_ELPA_INCLUDE_DIR}")
    message(STATUS "  ELPA_LIBRARIES=${ROKKO_ELPA_LIBRARIES}")
  endif(ROKKO_HAVE_ELPA)

  message(STATUS "  ROKKO_HAVE_ELEMENTAL:" ${ROKKO_HAVE_ELEMENTAL})
  if(ROKKO_HAVE_ELEMENTAL)
    message(STATUS "  ELEMENTAL_INCLUDE_DIR=${ROKKO_ELEMENTAL_INCLUDE_DIR}")
    message(STATUS "  ELEMENTAL_LIBRARIES=${ROKKO_ELEMENTAL_LIBRARIES}")
  endif(ROKKO_HAVE_ELEMENTAL)

  # parallel sparse solvers
  message(STATUS "  ROKKO_HAVE_PETSC:" ${ROKKO_HAVE_PETSC})
  if(ROKKO_HAVE_PETSC)
    message(STATUS "  PETSC_INCLUDE_DIR=${ROKKO_PETSC_INCLUDE_DIR}")
    message(STATUS "  PETSC_LIBRARIES=${ROKKO_PETSC_LIBRARIES}")
  endif(ROKKO_HAVE_PETSC)

  message(STATUS "  ROKKO_HAVE_SLEPC:" ${ROKKO_HAVE_SLEPC})
  if(ROKKO_HAVE_SLEPC)
    message(STATUS "  SLEPC_INCLUDE_DIR=${ROKKO_SLEPC_INCLUDE_DIR}")
    message(STATUS "  SLEPC_LIBRARIES=${ROKKO_SLEPC_LIBRARIES}")
  endif(ROKKO_HAVE_SLEPC)

  message(STATUS "  ROKKO_HAVE_ANASAZI:" ${ROKKO_HAVE_ANASAZI})
  if(ROKKO_HAVE_ANASAZI)
    message(STATUS "  ANASAZI_INCLUDE_DIR=${ROKKO_ANASAZI_INCLUDE_DIR}")
    message(STATUS "  ANASAZI_LIBRARIES=${ROKKO_ANASAZI_LIBRARIES}")
  endif(ROKKO_HAVE_ANASAZI)

  # paths for the all solvers
  message(STATUS "ROKKO_INCLUDE_DIR=${ROKKO_INCLUDE_DIR}")
  message(STATUS "ROKKO_LIBRARIES=${ROKKO_LIBRARIES}")

  message(STATUS "ROKKO_HAVE_SERIAL_DENSE_SOLVER:" ${ROKKO_HAVE_SERIAL_DENSE_SOLVER})
  message(STATUS "ROKKO_HAVE_PARALLEL_DENSE_SOLVER:" ${ROKKO_HAVE_PARALLEL_DENSE_SOLVER})
  message(STATUS "ROKKO_HAVE_PARALLEL_SPARSE_SOLVER:" ${ROKKO_HAVE_PARALLEL_SPARSE_SOLVER})

  # bindings
  message(STATUS "ROKKO_Python_BINDING:" ${ROKKO_Python_BINDING})

endmacro(display_rokko_components)
