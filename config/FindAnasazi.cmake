# - Try to find Anasazi
# Once done this will define
#
#  ANASAZI_FOUND        - system has Anasazi
#  ANASAZI_INCLUDE_DIR  - include directories for Anasazi
#  ANASAZI_LIBARIES     - libraries for Anasazi
#  ANASAZI_DIR          - directory where Anasazi is built

if(DEFINED ANASAZI_FOUND)
  return()
endif(DEFINED ANASAZI_FOUND)
  
message(STATUS "Checking for Anasazi library")
set(ANASAZI_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(ROKKO_SOLVER_DIR)
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR})
endif(ROKKO_SOLVER_DIR)
list(APPEND _PATHS ${CMAKE_INSTALL_PREFIX} "/opt/nano/rokko" "/opt/rokko" "/opt" "$ENV{HOME}/opt/rokko" "$ENV{HOME}/opt")

# Try to figure out ANASAZI_DIR by finding Anasazi_config.h
find_path(ANASAZI_DIR include/Anasazi_config.h
  PATHS ${_PATHS}
  DOC "Anasazi directory")

if(ANASAZI_DIR)
  set(ANASAZI_INCLUDE_DIR "${ANASAZI_DIR}/include")
else(ANASAZI_DIR)
  message(STATUS "Anasazi library: not found")
  return()
endif(ANASAZI_DIR)

set(_LIBS anasazitpetra ModeLaplace anasaziepetra anasazi thyraepetra thyracore tpetraext tpetrainout tpetra epetra kokkosdisttsqr kokkosnodetsqr kokkoslinalg kokkosnodeapi kokkos rtop tpi teuchosremainder teuchosnumerics teuchoscomm teuchosparameterlist teuchoscore)
foreach(name ${_LIBS})
  unset(_LIB CACHE)
  find_library(_LIB NAMES ${name} PATHS ${ANASAZI_DIR}/lib)
  if(_LIB)
    # message(STATUS ${name} ${_LIB})
    set(ANASAZI_FOUND TRUE)
    list(APPEND ANASAZI_LIBRARIES ${_LIB})
  endif(_LIB)
endforeach(name)
# mark_as_advanced(ANASAZI_LIBRARIES)

if(${ANASAZI_FOUND})
  message(STATUS "Anasazi include directory: ${ANASAZI_INCLUDE_DIR}")
  message(STATUS "Anasazi libraries: ${ANASAZI_LIBRARIES}")
else(${ANASAZI_FOUND})
  message(STATUS "Anasazi library: not found")
endif(${ANASAZI_FOUND})
