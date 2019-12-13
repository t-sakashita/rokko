#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(TESTS frank helmert heisenberg eigen_matrix eigen_matrix_fixed_size matrix_major matrix_major_template_transformation)
foreach(name ${TESTS})
  set(target_name test_generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${target_name} ${name})
endforeach(name)
unset(TESTS)

set(TESTS xyz)
foreach(name ${TESTS})
  set(target_name test_generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${target_name} ${name} ${CMAKE_SOURCE_DIR}/test/input_data/xyz_1triangle.ip)
endforeach(name)

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  set(TESTS heisenberg_mpi laplacian_mpi laplacian_mfree_mpi frank_mpi frank_function_pointer_mpi frank_functor_mpi eigen2distributed_matrix xyz_mpi)
  foreach(name ${TESTS})
    set(target_name test_generate_matrix_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
