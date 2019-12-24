#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# core wrappers
add_subdirectory(blas)
add_subdirectory(lapack)
if(ROKKO_HAVE_SCALAPACK)
  add_subdirectory(pblas)
  add_subdirectory(scalapack)
endif(ROKKO_HAVE_SCALAPACK)

add_subdirectory(mpi_unit_test)
add_subdirectory(solver_class)
add_subdirectory(grid)
add_subdirectory(mapping_bc)
add_subdirectory(generate_matrix)
add_subdirectory(diagonalize)
add_subdirectory(mpi_collective)
add_subdirectory(factory)
add_subdirectory(utility)

if(BUILD_C_BINDING)
  add_subdirectory(c)
endif(BUILD_C_BINDING)
