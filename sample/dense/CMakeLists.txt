#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

configure_file(xyz.dat.in ${CMAKE_CURRENT_BINARY_DIR}/xyz.dat COPYONLY)

set(SAMPLE dsyev)
foreach(name ${SAMPLE})
  set(target_name sample_dense_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${target_name} rokko)
  add_rokko_test(${name})
endforeach(name)
unset(SAMPLE)

set(SAMPLE frank heisenberg xyz)
foreach(name ${SAMPLE})
  set(target_name sample_dense_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  add_rokko_test(${name})
endforeach(name)
unset(SAMPLE)

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  set(SAMPLE frank_mpi frank_mpi_mapping heisenberg_mpi xyz_mpi all_gather_mpi gather_mpi)
  foreach(name ${SAMPLE})
    set(target_name sample_dense_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
  unset(SAMPLE)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
