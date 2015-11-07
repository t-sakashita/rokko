#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_subdirectory(dense_frank_mpi)

if(BUILD_PARALLEL_SPARSE_SOLVER)
  add_subdirectory(sparse_heisenberg_mpi)
endif(BUILD_PARALLEL_SPARSE_SOLVER)




