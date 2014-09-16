#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>,
#               2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# single process tests
set(TESTS localized_matrix matrix_major serial_dense_solvers solver)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)

# multi process tests
if(BUILD_ANASAZI)
  include_directories(${ANASAZI_INCLUDE_DIR})
endif(BUILD_ANASAZI)
set(TESTS grid localized2distributed_matrix distributed_matrix frank_functor_mpi frank_function_pointer_mpi mpi_thread_multiple parallel_dense_solvers parallel_sparse_solvers)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

set(TESTS grid_1d_mpi)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

# 2,3,4,5,6 processes

set(TESTS grid_mpi gather scatter)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-5 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

# add_rokko_test tests
set(TESTS timer)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
  add_rokko_test(${name})
endforeach(name)

# utility tests
set(TESTS sort_eigenpairs)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${Boost_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)

add_subdirectory(generate_matrix)
add_subdirectory(diagonalize)
