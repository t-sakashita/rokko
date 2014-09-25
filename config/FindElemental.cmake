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
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR} $ENV{ROKKO_SOLVER_DIR} ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE} ${CMAKE_INSTALL_PREFIX} $ENV{HOME}/opt/rokko/${CMAKE_BUILD_TYPE} $ENV{HOME}/opt/rokko $ENV{HOME}/opt/${CMAKE_BUILD_TYPE} $ENV{HOME}/opt /opt/rokko/${CMAKE_BUILD_TYPE} /opt/rokko /opt/${CMAKE_BUILD_TYPE} /opt)
endif(ELEMENTAL_DIR)

find_path(ELEMENTAL_DIR include/elemental.hpp PATHS ${_PATHS} DOC "Elemental directory")

if(ELEMENTAL_DIR)
  set(ELEMENTAL_INCLUDE_DIR "${ELEMENTAL_DIR}/include")
else(ELEMENTAL_DIR)
  message(STATUS "Elemental library: not found")
  set(ELEMENTAL_FOUND FALSE)
  return()
endif(ELEMENTAL_DIR)

find_library(_ELEMENTAL_LIBRARY
  NAMES elemental
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

find_library(_ELEMENTAL_LAPACK_ADDONS_LIBRARY
  NAMES lapack-addons
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal lapack_addons library")
if(_ELEMENTAL_LAPACK_ADDONS_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_LAPACK_ADDONS_LIBRARY})
endif(_ELEMENTAL_LAPACK_ADDONS_LIBRARY)

find_library(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY
  NAMES elem-dummy-lib
  PATHS ${ELEMENTAL_DIR}/lib
  DOC "The Elemetnal elem_dummy_lib library")
if(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY})
endif(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY)

set(ELEMENTAL_FOUND TRUE)
message(STATUS "Elemental include directory: ${ELEMENTAL_INCLUDE_DIR}")
message(STATUS "Elemental libraries: ${ELEMENTAL_LIBRARIES}")
