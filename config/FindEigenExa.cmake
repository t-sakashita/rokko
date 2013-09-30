# - Try to find EigenExa
#
# Once done this will define
#
#  EIGENEXA_FOUND        - system has EigenExa
#  EIGENEXA_INCLUDE_DIR  - include directories for EigenExa
#  EIGENEXA_LIBARIES     - libraries for EigenExa

if(DEFINED EIGENEXA_FOUND)
  return()
endif(DEFINED EIGENEXA_FOUND)
  
message(STATUS "Checking for EigenExa library")

# Standard search path
set(_PATHS "")
if(ROKKO_SOLVER_DIR)
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR})
endif(ROKKO_SOLVER_DIR)
list(APPEND _PATHS ${CMAKE_INSTALL_PREFIX} "/opt/nano/rokko" "/opt/rokko" "/opt" "$ENV{HOME}/opt/rokko" "$ENV{HOME}/opt")

find_path(_EIGENEXA_DIR include/eigen_libs.mod
  HINTS ${EIGENEXA_DIR} $ENV{EIGENEXA_DIR} PATHS ${_PATHS}
  DOC "EigenExa directory")
if(_EIGENEXA_DIR)
  set(EIGENEXA_INCLUDE_DIR "${_EIGENEXA_DIR}/include")
else(_EIGENEXA_DIR)
  message(STATUS "EigenExa library: not found")
  set(EIGENEXA_FOUND FALSE)
  return()
endif(_EIGENEXA_DIR)

find_library(_EIGENEXA_LIBRARY
  NAME EigenExa
  PATHS ${_EIGENEXA_DIR}/lib
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
