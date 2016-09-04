#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  add_subdirectory(dense_frank_mpi)
  add_subdirectory(dense_minij_mpi)
  add_subdirectory(dense_heisenberg_ladder_mpi)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)

if(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
  add_subdirectory(sparse_heisenberg_file_mpi)
  add_subdirectory(sparse_heisenberg_ladder_mpi)
endif(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)

