#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(TESTS frank)
if(BUILD_PARALLEL_DENSE_SOLVER)
  set(TESTS ${TESTS} frank_mpi frank_block_mpi)
endif(BUILD_PARALLEL_DENSE_SOLVER)
foreach(name ${TESTS})
  add_test(python_dense_${name}
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
