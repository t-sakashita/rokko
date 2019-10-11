#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(SAMPLE frank)
foreach(name ${SAMPLE})
  set(target_name example_c_dense_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${target_name} rokko)
endforeach(name)
unset(SAMPLE)

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  set(SAMPLE frank_mpi frank_block_mpi frank_function_pointer_mpi)
  foreach(name ${SAMPLE})
    set(target_name example_c_dense_${name})
    add_executable(${target_name} ${name}.c)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(${target_name} rokko)
  endforeach(name)
  unset(SAMPLE)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)


# gev_fixedB_mpi
if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  set(SAMPLE gev_fixedB_mpi)
  foreach(name ${SAMPLE})
    set(target_name example_c_dense_${name})
    add_executable(${target_name} ${name}.c gev_fixedB_mpi_c.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(${target_name} rokko)
  endforeach(name)
  unset(SAMPLE)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)

add_subdirectory(generate_matrix)
