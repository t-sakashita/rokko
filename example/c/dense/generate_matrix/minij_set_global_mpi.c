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
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_grid grid;
  struct rokko_mapping_bc map;
  struct rokko_distributed_matrix mat;
  int provided, myrank, nprocs;
  int local_i, local_j, global_i, global_j;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  dim = 10;
  if (myrank == 0) {
    printf("dimension = %d\n", dim);
  }

  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);
  rokko_mapping_bc_construct_block_size(&map, dim, 1, grid);
  rokko_distributed_matrix_construct(&mat, map);

  /* generate minij matrix */
  for(global_j=0; global_j<dim; ++global_j) {
    for(global_i=0; global_i<dim; ++global_i) {
      rokko_distributed_matrix_set_global(mat, global_i, global_j, MIN(global_i, global_j) + 1);
    }
  }
  rokko_distributed_matrix_print(mat);

  rokko_distributed_matrix_destruct(&mat);
  rokko_mapping_bc_destruct(&map);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
