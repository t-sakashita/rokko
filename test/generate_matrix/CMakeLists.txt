#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2021 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(TESTS frank minij laplacian helmert matrix012 heisenberg eigen_matrix eigen_matrix_fixed_size matrix_major matrix_major_template_transformation trace_mpi)
foreach(name ${TESTS})
  set(target_name test_generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${target_name} ${name})
endforeach(name)
unset(TESTS)

set(TESTS xyz xyz_file)
foreach(name ${TESTS})
  set(target_name test_generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${target_name} ${name} ${CMAKE_SOURCE_DIR}/test/input_data/xyz_1triangle.ip)
endforeach(name)

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  set(TESTS heisenberg_mpi xyz_mpi)
  foreach(name ${TESTS})
    set(target_name test_generate_matrix_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    foreach(np 1 2 4)
      add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
    endforeach(np)
  endforeach(name)
  set(TESTS laplacian_mpi laplacian_mfree_mpi frank_mpi frank_function_pointer_mpi frank_functor_mpi minij_mpi matrix012_mpi eigen2distributed_matrix rank_number)
  foreach(name ${TESTS})
    set(target_name test_generate_matrix_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko gtest_main)
    foreach(np 1 2 3 4 5 6)
      add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
    endforeach(np)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
