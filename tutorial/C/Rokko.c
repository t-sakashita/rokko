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

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_parallel_dense_ev solver;
  struct rokko_distributed_matrix frank, eigvecs;
  struct rokko_grid grid;
  struct rokko_mapping_bc map;
  struct rokko_eigen_vector eigvals;
  char* solver_name;

  int provided, ierr, myrank, nprocs, i, j;

  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  dim = 4;
  if (argc == 2) {
    solver_name = argv[1];
  } else {
    fprintf(stderr, "error: %s solver_name\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  printf("solver name = %s\n", solver_name);
  printf("matrix dimension = %d\n", dim);

  rokko_parallel_dense_ev_construct(&solver, solver_name, argc, argv);
  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);
  map = rokko_parallel_dense_ev_default_mapping(solver, dim, grid);
  rokko_distributed_matrix_construct(&frank, map);
  rokko_distributed_matrix_construct(&eigvecs, map);
  rokko_eigen_vector_construct(&eigvals, dim);

  /* generate frank matrix */
  for(i = 0; i<dim; ++i){
    for(j = 0; j<dim; ++j){
      double val = dim - (i>j?i:j);
      rokko_distributed_matrix_set_global(frank, i, j, val);
    }
  }
  rokko_distributed_matrix_print(frank);

  MPI_Barrier(MPI_COMM_WORLD);

  rokko_parallel_dense_ev_diagonalize_distributed_matrix(solver, frank, eigvals, eigvecs);

  MPI_Barrier(MPI_COMM_WORLD);

  if (myrank == 0) {
    printf("Eigenvalues: \n");
    for (i = 0; i < dim; ++i)
      printf("%30.20f\n", rokko_eigen_vector_get(eigvals, i));

    printf("Eigenstates: \n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  rokko_distributed_matrix_print(eigvecs);

  rokko_eigen_vector_destruct(&eigvals);
  rokko_distributed_matrix_destruct(&eigvecs);
  rokko_distributed_matrix_destruct(&frank);
  rokko_parallel_dense_ev_destruct(&solver);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
