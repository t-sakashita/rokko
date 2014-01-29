# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDE_DIR  - the PETSc include directory
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if(DEFINED PETSC_FOUND)
  return()
endif(DEFINED PETSC_FOUND)
  
message(STATUS "Checking for PETSc library")
set(PETSC_FOUND FALSE)

# Standard search path
set(_PATHS "")
if(ROKKO_SOLVER_DIR)
  list(APPEND _PATHS ${ROKKO_SOLVER_DIR})
endif(ROKKO_SOLVER_DIR)
list(APPEND _PATHS ${CMAKE_INSTALL_PREFIX} "/opt/nano/rokko" "/opt/rokko" "/opt" "$ENV{HOME}/opt/rokko" "$ENV{HOME}/opt")

# Standard paths for Debian with version number
file(GLOB tmp "/usr/lib/petscdir/*")
list(APPEND _PATHS ${tmp})
unset(tmp)

function (petsc_get_version)
  if (EXISTS "${PETSC_DIR}/include/petscversion.h")
    file (STRINGS "${PETSC_DIR}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" PARENT_SCOPE)
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" PARENT_SCOPE)
    endif ()
  else ()
    message (SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist")
  endif ()
endfunction ()

find_path(_PETSC_DIR include/petscversion.h
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR} PATHS ${_PATHS}
  DOC "PETSc directory")
if(_PETSC_DIR)
  set(PETSC_INCLUDE_DIR "${_PETSC_DIR}/include")
else(_PETSC_DIR)
  message(STATUS "Petsc library: not found")
  set(PETSC_FOUND FALSE)
  return()
endif(_PETSC_DIR)

find_library(_PETSC_LIBRARY
  NAMES petsc
  PATHS ${_PETSC_DIR}/lib
  DOC "The PETSC library")
if(_PETSC_LIBRARY)
  list(APPEND PETSC_LIBRARIES ${_PETSC_LIBRARY})
else(_PETSC_LIBRARY)
  message(STATUS "Petsc library: not found")
  set(PETSC_FOUND FALSE)
  return()
endif(_PETSC_LIBRARY)

set(PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc definitions" FORCE)

set(PETSC_FOUND TRUE)
message(STATUS "PETSc include directory: ${PETSC_INCLUDE_DIR}")
message(STATUS "PETSc libraries: ${PETSC_LIBRARIES}")
message(STATUS "PETSc definitions: ${PETSC_DEFINITIONS}")
