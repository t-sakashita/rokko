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

#pragma once

#include <stdbool.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  rokko_grid_col_major = 1, rokko_grid_row_major = 2
};

struct rokko_grid {
  void* ptr;
  int major;
};

void rokko_grid_construct(struct rokko_grid* grid, MPI_Comm comm, int grid_major);
void rokko_grid_construct_f(struct rokko_grid* grid, int comm, int grid_major);
void rokko_grid_destruct(struct rokko_grid* grid);
int rokko_grid_get_myrank(struct rokko_grid grid);
int rokko_grid_get_nprocs(struct rokko_grid grid);
int rokko_grid_get_myrow(struct rokko_grid grid);
int rokko_grid_get_mycol(struct rokko_grid grid);
int rokko_grid_get_nprow(struct rokko_grid grid);
int rokko_grid_get_npcol(struct rokko_grid grid);
bool rokko_grid_is_row_major(struct rokko_grid grid);
bool rokko_grid_is_col_major(struct rokko_grid grid);

#ifdef __cplusplus
}
#endif
