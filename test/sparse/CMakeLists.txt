#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2021 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
  set(TESTS laplacian_crs_mpi laplacian_mfree_mpi laplacian_crs_variable_local_rows_mpi laplacian_mfree_variable_local_rows_mpi tridiagonal_toeplitz_mfree_mpi distributed_crs_symmetric_matrix_extract)
  foreach(name ${TESTS})
    set(target_name test_sparse_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    foreach(np 1 2 3 4 5 6)
      add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
    endforeach(np)
  endforeach(name)
  set(TESTS heisenberg_mfree_output_market_through_crs)
  foreach(name ${TESTS})
    set(target_name test_sparse_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    foreach(np 1 2 4)
      add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
    endforeach(np)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
