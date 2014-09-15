/*****************************************************************************
 *
 * Rokko: Integrated Interface for libraries of eigenvalue decomposition
 *
 * Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
 *                            Synge Todo <wistaria@comp-phys.org>,
 *                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
 *    
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *****************************************************************************/

#include <mpi.h>
#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_parallel_dense_solver solver;
  struct rokko_distributed_matrix mat, Z;
  struct rokko_grid grid;
  struct rokko_localized_vector w;
  char* solver_name;

  int provided, ierr, myrank, nprocs, i;

  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (argc == 3) {
    solver_name = argv[1];
    dim = atoi(argv[2]);
  } else {
    fprintf(stderr, "error: %s solver_name dimension\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 34);
  }
    
  printf("solver name = %s\n", solver_name);
  printf("matrix dimension = %d\n", dim);

  rokko_parallel_dense_solver_construct(&solver, solver_name, argc, argv);
  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);

  rokko_distributed_matrix_construct(&mat, dim, dim, grid, solver, rokko_matrix_col_major);
  rokko_distributed_matrix_construct(&Z, dim, dim, grid, solver, rokko_matrix_col_major);
  rokko_localized_vector_construct(&w, dim);

  /* generate frank matrix */
  rokko_frank_matrix_generate_distributed_matrix(&mat);
  rokko_distributed_matrix_print(mat);

  rokko_parallel_dense_solver_diagonalize_distributed_matrix(&solver, &mat, &w, &Z);

  if (myrank == 0) {
    printf("Computed Eigenvalues =\n");
    for (i = 0; i < dim; ++i)
      printf("%30.20f\n", rokko_localized_vector_get(w, i));
  }

  rokko_distributed_matrix_destruct(&mat);
  rokko_distributed_matrix_destruct(&Z);
  rokko_localized_vector_destruct(&w);
  rokko_parallel_dense_solver_destruct(&solver);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
