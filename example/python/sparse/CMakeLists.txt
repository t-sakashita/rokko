#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(BUILD_PARALLEL_SPARSE_SOLVER)
  set(TESTS heisenberg laplacian)
endif(BUILD_PARALLEL_SPARSE_SOLVER)
foreach(name ${TESTS})
  pybind11_add_module(${name}_package py_${name}_mfree.cpp)
  target_include_directories(${name}_package PUBLIC ${MPI4PY_DIR}/include)
  target_link_libraries(${name}_package PRIVATE rokko)
  set_target_properties(${name}_package PROPERTIES OUTPUT_NAME ${name})
endforeach(name)


if(BUILD_PARALLEL_SPARSE_SOLVER)
  set(TESTS heisenberg_crs_mpi)
endif(BUILD_PARALLEL_SPARSE_SOLVER)
foreach(name ${TESTS})
  add_test(python_${name}
    ${CMAKE_COMMAND}
    -Dpython_interpreter=${PYTHON_EXECUTABLE}
    -Dcmd=${name}.py
    -Dinput=${name}
    -Doutput=${name}
    -Dpythonpath=${PROJECT_BINARY_DIR}/rokko:${PROJECT_SOURCE_DIR}/rokko
    -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
    -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
    -Dcmddir=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${PROJECT_SOURCE_DIR}/config/run_python_test.cmake
    )
endforeach(name)

