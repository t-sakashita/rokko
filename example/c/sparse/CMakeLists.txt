#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
  set(TESTS heisenberg_crs_mpi laplacian_crs_mpi heisenberg_mfree_mpi laplacian_mfree_mpi laplacian_mfree_no_context_mpi)
  foreach(name ${TESTS})
    set(target_name example_c_sparse_${name})
    add_executable(${target_name} ${name}.c)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(${target_name} rokko)
    add_test(${target_name}-1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${target_name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${target_name}-3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${target_name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${target_name}-5 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(${target_name}-6 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
