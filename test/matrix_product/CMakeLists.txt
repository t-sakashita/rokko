#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  set(TESTS dot_product_mpi product_mpi product_v_mpi)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} rokko gtest_main)
    add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-5 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endforeach(name)
  unset(TESTS)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
