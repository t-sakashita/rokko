# - Try to find ELPA
# Once done this will define
#
#  ELPA_FOUND        - system has ELPA
#  ELPA_LIBARIES     - libraries for ELPA

if(DEFINED ELPA_FOUND)
  return()
endif(DEFINED ELPA_FOUND)
  
message(STATUS "Checking for ELPA library")
set(ELPA_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(ELPA_DIR)
  set(_PATHS ${ELPA_DIR})
else(ELPA_DIR)
  list(APPEND _PATHS
  	      ${ELPA_ROOT}/${CMAKE_BUILD_TYPE}
	      ${ELPA_ROOT}
  	      $ENV{ELPA_ROOT}/${CMAKE_BUILD_TYPE}
	      $ENV{ELPA_ROOT}
  	      ${ROKKO_SOLVER_ROOT}/elpa/${CMAKE_BUILD_TYPE}
	      ${ROKKO_SOLVER_ROOT}/elpa
  	      $ENV{ROKKO_SOLVER_ROOT}/elpa/${CMAKE_BUILD_TYPE}
	      $ENV{ROKKO_SOLVER_ROOT}/elpa
	      ${CMAKE_INSTALL_PREFIX}/elpa/${CMAKE_BUILD_TYPE}
	      ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/elpa/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/elpa
	      /opt/rokko/elpa/${CMAKE_BUILD_TYPE}
	      /opt/rokko/elpa
	      /opt/rokko/${CMAKE_BUILD_TYPE}
	      /opt/rokko
	      /opt/local /opt
	      )
endif(ELPA_DIR)

foreach (_PATH ${_PATHS})
  list(APPEND _LIBPATHS "${_PATH}/lib")
endforeach()

find_library(_ELPA_LIBRARY
  NAMES elpa_openmp elpa
  PATHS ${_LIBPATHS}
  DOC "The ELPA library")
if(_ELPA_LIBRARY)
  list(APPEND ELPA_LIBRARIES ${_ELPA_LIBRARY})
else(_ELPA_LIBRARY)
  message(STATUS "ELPA library: not found")
  set(ELPA_FOUND FALSE)
  return()
endif(_ELPA_LIBRARY)

set(ELPA_FOUND TRUE)
message(STATUS "ELPA libraries: ${ELPA_LIBRARIES}")
