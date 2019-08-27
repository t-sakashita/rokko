# - Try to find mpi4py
# Once done this will define
#
#  MPI4PY_FOUND        - system has mpi4py
#  MPI4PY_DIR          - directory where mpi4py is installed

if(DEFINED MPI4PY_FOUND)
  return()
endif(DEFINED MPI4PY_FOUND)

message(STATUS "Checking for mpi4py for Python Binding")
set(MPI4PY_FOUND FALSE)

execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c"
"import mpi4py
print(mpi4py.__path__[0])
"
    RESULT_VARIABLE _PYTHON_SUCCESS
    OUTPUT_VARIABLE _PYTHON_OUTPUT_VALUES
    ERROR_VARIABLE _PYTHON_ERROR_VALUES
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

if(_PYTHON_SUCCESS MATCHES 0)
  set(MPI4PY_DIR ${_PYTHON_OUTPUT_VALUES})
  message(STATUS "Found mpi4py: ${MPI4PY_DIR}")
  set(MPI4PY_FOUND TRUE)
else(_PYTHON_SUCCESS MATCHES 0)
  message(FATAL_ERROR "Failed to find mpi4py:")
  message(FATAL_ERROR "  error: ${_PYTHON_ERROR_VALUES}")
  set(MPI4PY_FOUND FALSE)
endif(_PYTHON_SUCCESS MATCHES 0)
