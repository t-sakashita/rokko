#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# core wrappers
add_subdirectory(blas)
add_subdirectory(lapack)
if(ROKKO_HAVE_SCALAPACK)
  add_subdirectory(pblas)
  add_subdirectory(scalapack)
endif(ROKKO_HAVE_SCALAPACK)

# single process tests
set(TESTS serial_dense_ev solver)
foreach(name ${TESTS})
  set(target_name test_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)

# multi process tests
if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  if(BUILD_PARALLEL_DENSE_SOLVER)
    set(TESTS parallel_dense_ev)
  endif(BUILD_PARALLEL_DENSE_SOLVER)
  if(BUILD_PARALLEL_SPARSE_SOLVER)
    set(TESTS ${TESTS} parallel_sparse_ev)
  endif(BUILD_PARALLEL_SPARSE_SOLVER)
  foreach(name ${TESTS})
    set(target_name test_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    add_test(${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endforeach(name)
  unset(TESTS)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)


add_subdirectory(mpi_unit_test)
add_subdirectory(grid)
add_subdirectory(generate_matrix)
add_subdirectory(diagonalize)
add_subdirectory(mpi_collective)
add_subdirectory(factory)
add_subdirectory(utility)

if(BUILD_C_BINDING)
  add_subdirectory(c)
endif(BUILD_C_BINDING)
