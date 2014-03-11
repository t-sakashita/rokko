# - Try to find Elemental
#
# Once done this will define
#
#  ELEMENTAL_FOUND        - system has Elemental
#  ELEMENTAL_INCLUDE_DIR  - include directories for Elemental
#  ELEMENTAL_LIBARIES     - libraries for Elemental

if(DEFINED ELEMENTAL_FOUND)
  return()
endif(DEFINED ELEMENTAL_FOUND)
  
message(STATUS "Checking for Elemental library")

# Standard search path
set(_PATHS "")
if(ROKKO_SOLVER_DIR)
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR})
endif(ROKKO_SOLVER_DIR)
list(APPEND _PATHS ${CMAKE_INSTALL_PREFIX} "/opt/nano/rokko" "/opt/rokko" "/opt" "$ENV{HOME}/opt/rokko" "$ENV{HOME}/opt")

find_path(_ELEMENTAL_DIR include/elemental.hpp
  HINTS ${ELEMENTAL_DIR} $ENV{ELEMENTAL_DIR} PATHS ${_PATHS}
  DOC "Elemental directory")
if(_ELEMENTAL_DIR)
  set(ELEMENTAL_INCLUDE_DIR "${_ELEMENTAL_DIR}/include")
else(_ELEMENTAL_DIR)
  message(STATUS "Elemental library: not found")
  set(ELEMENTAL_FOUND FALSE)
  return()
endif(_ELEMENTAL_DIR)

find_library(_ELEMENTAL_LIBRARY
  NAMES elemental
  PATHS ${_ELEMENTAL_DIR}/lib
  DOC "The Elemetnal library")
if(_ELEMENTAL_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_LIBRARY})
else(_ELEMENTAL_LIBRARY)
  message(STATUS "Elemental library: not found")
  set(ELEMENTAL_FOUND FALSE)
  return()
endif(_ELEMENTAL_LIBRARY)

find_library(_ELEMENTAL_PMRRR_LIBRARY
  NAMES pmrrr
  PATHS ${_ELEMENTAL_DIR}/lib
  DOC "The Elemetnal pmrrr library")
if(_ELEMENTAL_PMRRR_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_PMRRR_LIBRARY})
endif(_ELEMENTAL_PMRRR_LIBRARY)

find_library(_ELEMENTAL_LAPACK_ADDONS_LIBRARY
  NAMES lapack-addons
  PATHS ${_ELEMENTAL_DIR}/lib
  DOC "The Elemetnal lapack_addons library")
if(_ELEMENTAL_LAPACK_ADDONS_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_LAPACK_ADDONS_LIBRARY})
endif(_ELEMENTAL_LAPACK_ADDONS_LIBRARY)

find_library(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY
  NAMES elem-dummy-lib
  PATHS ${_ELEMENTAL_DIR}/lib
  DOC "The Elemetnal elem_dummy_lib library")
if(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY)
  list(APPEND ELEMENTAL_LIBRARIES ${_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY})
endif(_ELEMENTAL_ELEM_DUMMY_LIB_LIBRARY)

set(ELEMENTAL_FOUND TRUE)
message(STATUS "Elemental include directory: ${ELEMENTAL_INCLUDE_DIR}")
message(STATUS "Elemental libraries: ${ELEMENTAL_LIBRARIES}")
