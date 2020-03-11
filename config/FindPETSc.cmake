# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDE_DIR  - the PETSc include directory
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_DIR          - directory where PETSc is installed

if(DEFINED PETSC_FOUND)
  return()
endif(DEFINED PETSC_FOUND)
  
message(STATUS "Checking for PETSc library")
set(PETSC_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(PETSC_DIR)
  set(_PATHS ${PETSC_DIR})
else(PETSC_DIR)
  list(APPEND _PATHS
  	      ${PETSC_ROOT}/${CMAKE_BUILD_TYPE}
	      ${PETSC_ROOT}
  	      $ENV{PETSC_ROOT}/${CMAKE_BUILD_TYPE}
	      $ENV{PETSC_ROOT}
  	      ${ROKKO_SOLVER_ROOT}/petsc/${CMAKE_BUILD_TYPE}
	      ${ROKKO_SOLVER_ROOT}/petsc
  	      $ENV{ROKKO_SOLVER_ROOT}/petsc/${CMAKE_BUILD_TYPE}
	      $ENV{ROKKO_SOLVER_ROOT}/petsc
	      ${CMAKE_INSTALL_PREFIX}/petsc/${CMAKE_BUILD_TYPE}
	      ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/petsc/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/petsc
	      /opt/rokko/petsc/${CMAKE_BUILD_TYPE}
	      /opt/rokko/petsc
	      /opt/rokko/${CMAKE_BUILD_TYPE}
	      /opt/rokko
	      /opt/local /opt
              /usr
	      )
  # Standard paths for Debian with version number
  file(GLOB tmp "/usr/lib/petscdir/*")
  list(APPEND _PATHS ${tmp})
  unset(tmp)
endif(PETSC_DIR)

find_path(PETSC_DIR NAMES include/petscversion.h PATHS ${_PATHS} DOC "PETSc directory")
if(PETSC_DIR)
  set(PETSC_INCLUDE_DIR "${PETSC_DIR}/include")
else(PETSC_DIR)
  find_path(PETSC_DIR NAMES include/petsc/petscversion.h PATHS ${_PATHS} DOC "PETSc directory")
  if(PETSC_DIR)
    set(PETSC_INCLUDE_DIR "${PETSC_DIR}/include/petsc")
  else(PETSC_DIR)
    message(STATUS "PETSc library: include file not found")
    set(PETSC_FOUND FALSE)
    return()
  endif(PETSC_DIR)
endif(PETSC_DIR)
message(STATUS "PETSc include directory: ${PETSC_INCLUDE_DIR}")

find_library(_PETSC_LIBRARY
  NAMES petsc
  PATHS ${PETSC_DIR}/lib
  DOC "The PETSC library")
if(_PETSC_LIBRARY)
  list(APPEND PETSC_LIBRARIES ${_PETSC_LIBRARY})
else(_PETSC_LIBRARY)
  message(STATUS "PETSc library: library file not found")
  set(PETSC_FOUND FALSE)
  return()
endif(_PETSC_LIBRARY)

set(PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc definitions" FORCE)

set(PETSC_FOUND TRUE)
message(STATUS "PETSc libraries: ${PETSC_LIBRARIES}")
message(STATUS "PETSc definitions: ${PETSC_DEFINITIONS}")
