# - Try to find Elpa
# Once done this will define
#
#  ELPA_FOUND        - system has Elpa
#  ELPA_LIBARIES     - libraries for Elpa

if(DEFINED ELPA_FOUND)
  return()
endif(DEFINED ELPA_FOUND)
  
message(STATUS "Checking for Elpa library")
set(ELPA_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(ELPA_DIR)
  set(_PATHS ${ELPA_DIR})
else(ELPA_DIR)
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR} $ENV{ROKKO_SOLVER_DIR} ${CMAKE_INSTALL_PREFIX} "$ENV{HOME}/opt/rokko" "$ENV{HOME}/opt" "/opt/rokko" "/opt")
endif(ELPA_DIR)

foreach (_PATH ${_PATHS})
  list(APPEND _LIBPATHS "${_PATH}/lib")
endforeach()

find_library(_ELPA_LIBRARY
  NAME elpa-2011.12
  PATHS ${_LIBPATHS}
  DOC "The Elpa library")
if(_ELPA_LIBRARY)
  list(APPEND ELPA_LIBRARIES ${_ELPA_LIBRARY})
else(_ELPA_LIBRARY)
  message(STATUS "Elpa library: not found")
  set(ELPA_FOUND FALSE)
  return()
endif(_ELPA_LIBRARY)

set(ELPA_FOUND TRUE)
message(STATUS "Elpa libraries: ${ELPA_LIBRARIES}")
