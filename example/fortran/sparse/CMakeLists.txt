#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)
include_directories(${MPI_Fortran_INCLUDE_PATH})

if(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
  set(SAMPLE heisenberg_crs_mpi laplacian_mfree_mpi laplacian_crs_mpi)
  foreach(name ${SAMPLE})
    set(target_name example_fortran_sparse_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko ${MPI_Fortran_LIBRARIES})
    if(ROKKO_USE_CXX_LINKER)
      SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    endif(ROKKO_USE_CXX_LINKER)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)

