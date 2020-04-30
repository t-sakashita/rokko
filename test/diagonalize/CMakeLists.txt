#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# copy files
set(PATTERN_FILES "heisenberg.ip")
foreach( pattern_file ${PATTERN_FILES} )
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file} ${CMAKE_CURRENT_BINARY_DIR}/${pattern_file} COPYONLY)
endforeach( pattern_file )

set(TESTS frank minij laplacian tridiagonal_toeplitz helmert heisenberg xyz xyz_file)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko gtest_main)
  add_rokko_test(${name})
endforeach(name)
unset(TESTS)

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  if(BUILD_PARALLEL_DENSE_SOLVER)
    set(TESTS frank_mpi minij_mpi laplacian_mpi tridiagonal_toeplitz_mpi helmert_mpi)
    foreach(name ${TESTS})
      add_executable(${name} ${name}.cpp)
      target_link_libraries(${name} rokko gtest_main)
      add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
      add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    endforeach(name)
    unset(TESTS)
  endif(BUILD_PARALLEL_DENSE_SOLVER)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
