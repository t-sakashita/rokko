#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  set(SAMPLE minij_mpi minij_block_size_mpi minij_mpi_thread_single)
  foreach(name ${SAMPLE})
    set(target_name use_rokko_dense_minij_mpi_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
  unset(SAMPLE)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)

