#  Copyright Olivier Parcollet and Matthias Troyer 2010.
#  Copyright Synge Todo 2015.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
#  Python settings : 
#
#  This module checks that : 
#  - the python interpreter is working and version >= 2.6
#  - it has modules : distutils, numpy, tables, scipy
# 
#  This module defines the variables
#  - PYTHON_FOUND : Was Python found
#  - PYTHON_EXECUTABLE : name of the python interpreter
#  - PYTHON_INCLUDE_DIRS : include for compilation
#  - PYTHON_LIBRARIES : link flags 

if (NOT PYTHON_EXECUTABLE)
  find_program(PYTHON_EXECUTABLE python PATHS $ENV{PATH})
  if (NOT PYTHON_EXECUTABLE)
    set (PYTHON_FOUND FALSE)
  else(NOT PYTHON_EXECUTABLE)
    set(PYTHON_FOUND TRUE)
  endif(NOT PYTHON_EXECUTABLE)
else (NOT PYTHON_EXECUTABLE)
  set(PYTHON_FOUND TRUE)
endif (NOT PYTHON_EXECUTABLE)

set(PYTHON_MINIMAL_VERSION 2.7)

if (WIN32)
  MESSAGE (STATUS "Looking for PythonLibs")
  find_package(PythonLibs)
endif (WIN32)

IF (PYTHON_FOUND)
  MESSAGE(STATUS "Python interpreter: ${PYTHON_EXECUTABLE}")
  FUNCTION(EXEC_PYTHON_SCRIPT the_script output_var_name)
    IF ("${PYTHON_EXECUTABLE}" MATCHES ".*ipython.*")
      EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} "--c=${the_script}" 
        OUTPUT_VARIABLE re    s RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
    ELSE ("${PYTHON_EXECUTABLE}" MATCHES ".*ipython.*")
      EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "${the_script}" 
        OUTPUT_VARIABLE res RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
    ENDIF ("${PYTHON_EXECUTABLE}" MATCHES ".*ipython.*")
    IF (NOT returncode EQUAL 0)
      MESSAGE(FATAL_ERROR "The script : ${the_script} \n did not run properly in the Python interpreter. Check your python installation.") 
    ENDIF (NOT returncode EQUAL 0)
    SET( ${output_var_name} ${res} PARENT_SCOPE)
  ENDFUNCTION (EXEC_PYTHON_SCRIPT)

  #
  # Check the interpreter and its version
  #
  EXEC_PYTHON_SCRIPT ("import sys, string; print sys.version.split()[0]" PYTHON_VERSION)
  STRING(COMPARE GREATER ${PYTHON_MINIMAL_VERSION} ${PYTHON_VERSION} PYTHON_VERSION_NOT_OK)
  IF (PYTHON_VERSION_NOT_OK)
    MESSAGE(WARNING "Python intepreter version is ${PYTHON_VERSION} . It should be >= ${PYTHON_MINIMAL_VERSION}")
    SET(PYTHON_FOUND FALSE)
  ENDIF (PYTHON_VERSION_NOT_OK)
  MESSAGE(STATUS "Python version: ${PYTHON_VERSION}")

  #
  # Check for Python include path
  #
  EXEC_PYTHON_SCRIPT ("import distutils ; from distutils.sysconfig import * ; print distutils.sysconfig.get_python_inc()"  PYTHON_INCLUDE_DIRS )
  message(STATUS "PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}" )
  mark_as_advanced(PYTHON_INCLUDE_DIRS)
  FIND_PATH(TEST_PYTHON_INCLUDE patchlevel.h PATHS ${PYTHON_INCLUDE_DIRS} NO_DEFAULT_PATH)
  if (NOT TEST_PYTHON_INCLUDE)
    message (ERROR "The Python header files have not been found. Please check that you installed the Python headers and not only the interpreter.")
  endif (NOT TEST_PYTHON_INCLUDE)

  #
  # Check for Python library path
  #
  EXEC_PYTHON_SCRIPT ("import string; from distutils.sysconfig import *; print '%s/config' % get_python_lib(0,1)" PYTHON_LIBRARY_BASE_PATH)
  EXEC_PYTHON_SCRIPT ("import string; from distutils.sysconfig import *; print 'libpython%s' % string.join(get_config_vars('VERSION'))" PYTHON_LIBRARY_BASE_FILE)
  IF(BUILD_SHARED_LIBS)
    FIND_FILE(PYTHON_LIBRARIES NAMES "${PYTHON_LIBRARY_BASE_FILE}.so" "${PYTHON_LIBRARY_BASE_FILE}.dylib" PATHS ${PYTHON_LIBRARY_BASE_PATH})
    IF(NOT PYTHON_LIBRARIES)
      FIND_FILE(PYTHON_LIBRARIES NAMES "${PYTHON_LIBRARY_BASE_FILE}.a" PATHS ${PYTHON_LIBRARY_BASE_PATH})
    ENDIF(NOT PYTHON_LIBRARIES)
  ELSE(BUILD_SHARED_LIBS)
    FIND_FILE(PYTHON_LIBRARIES NAMES "${PYTHON_LIBRARY_BASE_FILE}.a" PATHS ${PYTHON_LIBRARY_BASE_PATH})
  ENDIF(BUILD_SHARED_LIBS)
  MESSAGE(STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}" )
  mark_as_advanced(PYTHON_LIBRARIES)
ENDIF (PYTHON_FOUND)

# PYTHON_ADD_MODULE(<name> src1 src2 ... srcN) is used to build modules for python.
# PYTHON_WRITE_MODULES_HEADER(<filename>) writes a header file you can include
# in your sources to initialize the static python modules
function(PYTHON_ADD_MODULE _NAME )
  get_property(_TARGET_SUPPORTS_SHARED_LIBS
    GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
  option(PYTHON_ENABLE_MODULE_${_NAME} "Add module ${_NAME}" TRUE)
  option(PYTHON_MODULE_${_NAME}_BUILD_SHARED
    "Add module ${_NAME} shared" ${_TARGET_SUPPORTS_SHARED_LIBS})

  # Mark these options as advanced
  mark_as_advanced(PYTHON_ENABLE_MODULE_${_NAME}
    PYTHON_MODULE_${_NAME}_BUILD_SHARED)

  if(PYTHON_ENABLE_MODULE_${_NAME})
    if(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      set(PY_MODULE_TYPE MODULE)
    else()
      set(PY_MODULE_TYPE STATIC)
      set_property(GLOBAL  APPEND  PROPERTY  PY_STATIC_MODULES_LIST ${_NAME})
    endif()

    set_property(GLOBAL  APPEND  PROPERTY  PY_MODULES_LIST ${_NAME})
    add_library(${_NAME} ${PY_MODULE_TYPE} ${ARGN})
#    target_link_libraries(${_NAME} ${PYTHON_LIBRARIES})

    if(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      set_target_properties(${_NAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
      if(WIN32 AND NOT CYGWIN)
        set_target_properties(${_NAME} PROPERTIES SUFFIX ".pyd")
      endif()
    endif()

  endif()
endfunction()
