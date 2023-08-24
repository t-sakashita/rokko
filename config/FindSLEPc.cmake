# - Try to find SLEPC
# Once done this will define
#
#  SLEPC_FOUND        - system has SLEPc
#  SLEPC_INCLUDE_DIR  - include directories for SLEPc
#  SLEPC_LIBARIES     - libraries for SLEPc
#  SLEPC_DIR          - directory where SLEPc is installed
#
# Assumes that PETSC_DIR and PETSC_ARCH has been set by
# alredy calling find_package(PETSc)

if(DEFINED SLEPC_FOUND)
  return()
endif(DEFINED SLEPC_FOUND)

message(STATUS "Checking for SLEPc library")
find_package(PETSc)
set(SLEPC_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(SLEPC_DIR)
  set(_PATHS ${SLEPC_DIR})
else(SLEPC_DIR)
  list(APPEND _PATHS
  	      ${SLEPC_ROOT}/${CMAKE_BUILD_TYPE}
	      ${SLEPC_ROOT}
  	      $ENV{SLEPC_ROOT}/${CMAKE_BUILD_TYPE}
	      $ENV{SLEPC_ROOT}
  	      ${ROKKO_SOLVER_ROOT}/slepc/${CMAKE_BUILD_TYPE}
	      ${ROKKO_SOLVER_ROOT}/slepc
  	      $ENV{ROKKO_SOLVER_ROOT}/slepc/${CMAKE_BUILD_TYPE}
	      $ENV{ROKKO_SOLVER_ROOT}/slepc
	      ${CMAKE_INSTALL_PREFIX}/slepc/${CMAKE_BUILD_TYPE}
	      ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/slepc/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/slepc
	      /opt/rokko/slepc/${CMAKE_BUILD_TYPE}
	      /opt/rokko/slepc
	      /opt/rokko/${CMAKE_BUILD_TYPE}
	      /opt/rokko
	      /opt/local /opt
              /usr
	      )
  # Standard paths for Debian with version number
  file(GLOB tmp "/usr/lib/slepcdir/*")
  list(APPEND _PATHS ${tmp})
  unset(tmp)
endif(SLEPC_DIR)

# Try to figure out SLEPC_DIR by finding slepc.h
find_path(SLEPC_DIR NAMES include/slepc.h PATHS ${_PATHS} DOC "SLEPc directory")

if(SLEPC_DIR)
  set(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include")
else(SLEPC_DIR)
  find_path(SLEPC_DIR NAMES include/slepc/slepc.h PATHS ${_PATHS} DOC "SLEPc directory")
  if(SLEPC_DIR)
    set(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include/slepc")
  else(SLEPC_DIR)
    message(STATUS "SLEPc library: include file not found")
    set(SLEPC_FOUND FALSE)
    return()
  endif(SLEPC_DIR)
endif(SLEPC_DIR)
message(STATUS "SLEPc include directory: ${SLEPC_INCLUDE_DIR}")

find_library(SLEPC_LIBRARIES
  NAMES slepc
  HINTS ${SLEPC_DIR}/lib
  HINTS ${SLEPC_DIR}/${PETSC_ARCH}/lib
  DOC "The SLEPc library"
  )
mark_as_advanced(SLEPC_LIBRARIES)

set(SLEPC_FOUND TRUE)
message(STATUS "SLEPc libraries: ${SLEPC_LIBRARIES}")

return()
