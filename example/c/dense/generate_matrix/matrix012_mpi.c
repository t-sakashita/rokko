/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <rokko/rokko.h>
#include <rokko/utility/matrix012.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_grid grid;
  struct rokko_mapping_bc map;
  struct rokko_distributed_matrix mat;

  int provided, myrank, nprocs;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (myrank == 0) {
    printf("dimension = %d\n", dim);
  }

  dim = 10;
  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);
  rokko_mapping_bc_construct_block_size(&map, dim, 1, grid);
  rokko_distributed_matrix_construct(&mat, map);

  /* generate matrix012 */
  rokko_matrix012_generate_distributed_matrix(mat);
  rokko_distributed_matrix_print(mat);

  rokko_distributed_matrix_destruct(&mat);
  rokko_mapping_bc_destruct(&map);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
