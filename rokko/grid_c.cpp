/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/grid.hpp>
#include <rokko/rokko_dense.h>

void rokko_grid_construct(rokko_grid* grid, MPI_Comm comm, int grid_major) {
  if (grid_major == rokko_grid_col_major)
    grid->ptr = static_cast<void*>(new rokko::grid(comm, rokko::grid_col_major));
  else
    grid->ptr = static_cast<void*>(new rokko::grid(comm, rokko::grid_row_major));
  grid->major = grid_major;
}

void rokko_grid_construct_f(rokko_grid* grid, int comm_f, int grid_major) {
  MPI_Comm comm = MPI_Comm_f2c(comm_f);
  rokko_grid_construct(grid, comm, grid_major);
}

void rokko_grid_destruct(rokko_grid* grid) {
  delete static_cast<rokko::grid*>(grid->ptr);
  grid->ptr = 0;
}

int rokko_grid_get_myrank(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_myrank();
}

int rokko_grid_get_nprocs(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_nprocs();
}
