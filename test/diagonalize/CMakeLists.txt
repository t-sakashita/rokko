#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2021 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# copy files
set(PATTERN_FILES "heisenberg.ip")
foreach(pattern_file ${PATTERN_FILES})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file} ${CMAKE_CURRENT_BINARY_DIR}/${pattern_file} COPYONLY)
endforeach( pattern_file )

set(TESTS frank minij laplacian tridiagonal_toeplitz helmert heisenberg xyz xyz_file)
foreach(name ${TESTS})
  set(target_name test_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko gtest_main)
  add_test(${target_name} ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  if(BUILD_PARALLEL_DENSE_SOLVER)
    set(TESTS frank_mpi minij_mpi laplacian_mpi tridiagonal_toeplitz_mpi helmert_mpi)
    foreach(name ${TESTS})
      set(target_name test_diagonalize_${name})
      add_executable(${target_name} ${name}.cpp)
      set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
      target_link_libraries(${target_name} rokko gtest_main)
      foreach(np 1 2 3 4 5 6)
        add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
      endforeach(np)
    endforeach(name)
    unset(TESTS)
  endif(BUILD_PARALLEL_DENSE_SOLVER)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
