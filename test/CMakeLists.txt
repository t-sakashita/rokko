#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>,
#               2013-2013    Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Tests using Boost.test
# Check and build Boost.Test Unit Test Framework
find_package(Boost 1.48.0 COMPONENTS unit_test_framework)
if (Boost_FOUND)
  add_definitions(-DBOOST_TEST_DYN_LINK)
endif (Boost_FOUND)

# single process tests
set(TESTS solver_factory localized_matrix serial_diagonalize)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)

# multi process tests
set(TESTS grid solver distributed_matrix)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

set(TESTS grid_mpi)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

# 2 & 4 processes

set(TESTS diagonalize gather scatter product spin_hamiltonian_parallel_test)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

if (SLEPC_LIBRARIES)
  include_directories(${SLEPC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})
  set(TESTS xyz)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} rokko ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
    add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endforeach(name)
  unset(TESTS)
endif(SLEPC_LIBRARIES)

# add_rokko_test tests
set(TESTS frank timer spin_hamiltonian spin_hamiltonian2)
# set(TESTS frank scatter gather)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
  add_rokko_test(${name})
endforeach(name)

if(BUILD_LAPACKE)
  set(TESTS dsyev)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.c)
    target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
    add_rokko_test(${name})
  endforeach(name)
endif(BUILD_LAPACKE)
