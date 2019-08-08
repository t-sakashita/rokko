/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/grid.hpp>
#include <rokko/grid.h>

#ifdef __cplusplus
extern "C" {
#endif

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
  grid->ptr = nullptr;
}

int rokko_grid_get_myrank(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_myrank();
}

int rokko_grid_get_nprocs(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_nprocs();
}

int rokko_grid_get_myrow(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_myrow();
}

int rokko_grid_get_mycol(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_mycol();
}

int rokko_grid_get_nprow(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_nprow();
}

int rokko_grid_get_npcol(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->get_npcol();
}

bool rokko_grid_is_row_major(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->is_row_major();
}

bool rokko_grid_is_col_major(rokko_grid grid) {
  return static_cast<rokko::grid*>(grid.ptr)->is_col_major();
}

#ifdef __cplusplus
}
#endif

