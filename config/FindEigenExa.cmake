# - Try to find EigenExa
# Once done this will define
#
#  EIGENEXA_FOUND        - system has EigenExa
#  EIGENEXA_INCLUDE_DIR  - include directories for EigenExa
#  EIGENEXA_LIBARIES     - libraries for EigenExa
#  EIGENEXA_DIR          - directory where EigenExa is installed

if(DEFINED EIGENEXA_FOUND)
  return()
endif(DEFINED EIGENEXA_FOUND)
  
message(STATUS "Checking for EigenExa library")
set(EIGENEXA_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(EIGENEXA_DIR)
  set(_PATHS ${EIGENEXA_DIR})
else(EIGENEXA_DIR)
  list(APPEND _PATHS
  	      ${EIGENEXA_ROOT}/${CMAKE_BUILD_TYPE}
	      ${EIGENEXA_ROOT}
  	      $ENV{EIGENEXA_ROOT}/${CMAKE_BUILD_TYPE}
	      $ENV{EIGENEXA_ROOT}
  	      ${ROKKO_SOLVER_ROOT}/eigenexa/${CMAKE_BUILD_TYPE}
	      ${ROKKO_SOLVER_ROOT}/eigenexa
  	      $ENV{ROKKO_SOLVER_ROOT}/eigenexa/${CMAKE_BUILD_TYPE}
	      $ENV{ROKKO_SOLVER_ROOT}/eigenexa
	      ${CMAKE_INSTALL_PREFIX}/eigenexa/${CMAKE_BUILD_TYPE}
	      ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/eigenexa/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/eigenexa
	      /opt/rokko/eigenexa/${CMAKE_BUILD_TYPE}
	      /opt/rokko/eigenexa
	      /opt/rokko/${CMAKE_BUILD_TYPE}
	      /opt/rokko
	      /opt/local /opt
	      )
endif(EIGENEXA_DIR)

find_path(EIGENEXA_DIR include/eigen_libs_mod.mod PATHS ${_PATHS} DOC "EigenExa directory")

if(EIGENEXA_DIR)
  set(EIGENEXA_INCLUDE_DIR "${EIGENEXA_DIR}/include")
else(EIGENEXA_DIR)
  message(STATUS "EigenExa library: not found")
  return()
endif(EIGENEXA_DIR)

find_library(_EIGENEXA_LIBRARY
  NAME EigenExa
  PATHS ${EIGENEXA_DIR}/lib
  DOC "The EigenExa library")
if(_EIGENEXA_LIBRARY)
  list(APPEND EIGENEXA_LIBRARIES ${_EIGENEXA_LIBRARY})
else(_EIGENEXA_LIBRARY)
  message(STATUS "EigenExa library: not found")
  set(EIGENEXA_FOUND FALSE)
  return()
endif(_EIGENEXA_LIBRARY)

set(EIGENEXA_FOUND TRUE)
message(STATUS "EigenExa include directory: ${EIGENEXA_INCLUDE_DIR}")
message(STATUS "EigenExa libraries: ${EIGENEXA_LIBRARIES}")
