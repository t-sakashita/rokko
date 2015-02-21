#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# copy files
set(PATTERN_FILES "heisenberg.ip" "xyz_slepc_matfree.ip")
foreach( pattern_file ${PATTERN_FILES} )
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file} ${CMAKE_CURRENT_BINARY_DIR}/${pattern_file} COPYONLY)
endforeach( pattern_file )

set(TESTS frank heisenberg xyz xyz_file)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES})
  add_rokko_test(${name})
endforeach(name)
unset(TESTS)

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  set(TESTS dot_product_mpi product_mpi product_v_mpi)
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

  if(BUILD_PARALLEL_DENSE_SOLVER)
    set(TESTS frank_mpi)
    foreach(name ${TESTS})
      add_executable(${name} ${name}.cpp)
      target_link_libraries(${name} rokko ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
      add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
      add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    endforeach(name)
    unset(TESTS)

    if(BUILD_SLEPC)
      set(TESTS xyz_slepc_matfree)
        foreach(name ${TESTS})
        add_executable(${name} ${name}.cpp)
        target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES})
        add_rokko_test(${name})
      endforeach(name)
      unset(TESTS)
    endif(BUILD_SLEPC)
  endif(BUILD_PARALLEL_DENSE_SOLVER)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
