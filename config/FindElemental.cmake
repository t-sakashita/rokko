# - Try to find Elemental
# Once done this will define
#
#  ELEMENTAL_FOUND        - system has Elemental
#  ELEMENTAL_INCLUDE_DIR  - include directories for Elemental
#  ELEMENTAL_LIBARIES     - libraries for Elemental
#  ELEMENTAL_DIR          - directory where Elemental is installed

if(DEFINED ELEMENTAL_FOUND)
  return()
endif(DEFINED ELEMENTAL_FOUND)
  
message(STATUS "Checking for Elemental library")
set(ELEMENTAL_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(ELEMENTAL_DIR)
  set(_PATHS ${ELEMENTAL_DIR})
else(ELEMENTAL_DIR)
  list(APPEND _PATHS
  	      ${ELEMENTAL_ROOT}/${CMAKE_BUILD_TYPE}
	      ${ELEMENTAL_ROOT}
  	      $ENV{ELEMENTAL_ROOT}/${CMAKE_BUILD_TYPE}
	      $ENV{ELEMENTAL_ROOT}
  	      ${ROKKO_SOLVER_ROOT}/elemental/${CMAKE_BUILD_TYPE}
	      ${ROKKO_SOLVER_ROOT}/elemental
  	      $ENV{ROKKO_SOLVER_ROOT}/elemental/${CMAKE_BUILD_TYPE}
	      $ENV{ROKKO_SOLVER_ROOT}/elemental
	      ${CMAKE_INSTALL_PREFIX}/elemental/${CMAKE_BUILD_TYPE}
	      ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/elemental/${CMAKE_BUILD_TYPE}
	      $ENV{HOME}/rokko/elemental
	      /opt/rokko/elemental/${CMAKE_BUILD_TYPE}
	      /opt/rokko/elemental
	      /opt/rokko/${CMAKE_BUILD_TYPE}
	      /opt/rokko
	      /opt/local /opt
	      )
endif(ELEMENTAL_DIR)

find_path(ELEMENTAL_DIR include/El.h include/elemental.hpp PATHS ${_PATHS} DOC "Elemental directory")

if(ELEMENTAL_DIR)
  set(ELEMENTAL_INCLUDE_DIR "${ELEMENTAL_DIR}/include")
else(ELEMENTAL_DIR)
  message(STATUS "Elemental library: not found")
  set(ELEMENTAL_FOUND FALSE)
  return()
endif(ELEMENTAL_DIR)

find_library(_ELEMENTAL_LIBRARY
  NAMES El elemental
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal library")
if(_ELEMENTAL_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_LIBRARY})
else(_ELEMENTAL_LIBRARY)
  message(STATUS "Elemental library: not found")
  return()
endif(_ELEMENTAL_LIBRARY)

find_library(_ELEMENTAL_PMRRR_LIBRARY
  NAMES pmrrr
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal pmrrr library")
if(_ELEMENTAL_PMRRR_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_PMRRR_LIBRARY})
endif(_ELEMENTAL_PMRRR_LIBRARY)

find_library(_ELEMENTAL_KISS_FFT_LIBRARY
  NAMES kiss_fft
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal kiss_fft library")
if(_ELEMENTAL_KISS_FFT_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_KISS_FFT_LIBRARY})
endif(_ELEMENTAL_KISS_FFT_LIBRARY)

find_library(_ELEMENTAL_METIS_LIBRARY
  NAMES metis
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal metis library")
if(_ELEMENTAL_METIS_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_METIS_LIBRARY})
endif(_ELEMENTAL_METIS_LIBRARY)

set(ELEMENTAL_FOUND TRUE)
message(STATUS "Elemental include directory: ${ELEMENTAL_INCLUDE_DIR}")
message(STATUS "Elemental libraries: ${ELEMENTAL_LIBRARIES}")
