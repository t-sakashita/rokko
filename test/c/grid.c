/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  int i, n;
  int provided, ierr;
  struct rokko_grid grid;
  bool flag;
  int flag_int;
  
  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);

  flag = rokko_grid_is_row_major(grid);
  printf("rokko_grid_is_row_major = %d\n", flag);
  flag = rokko_grid_is_col_major(grid);
  printf("rokko_grid_is_col_major = %d\n", flag);

  flag_int = rokko_grid_is_row_major(grid);
  printf("rokko_grid_is_row_major = %d\n", flag_int);
  flag_int = rokko_grid_is_col_major(grid);
  printf("rokko_grid_is_col_major = %d\n", flag_int);

  rokko_grid_destruct(&grid);

  MPI_Finalize();
}
