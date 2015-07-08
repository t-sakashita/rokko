/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
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
  int provided, ierr, myrank, nprocs;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int root = 0;
  char default_solver_name[] = "anasazi";
  char* solver_name;
  if (argc >= 2) {
    solver_name = argv[1];
  } else {
    solver_name = default_solver_name;
  }
  int dim;
  if (argc == 3) {
    dim = atoi(argv[2]);
  } else {
    dim = 10;
  }

  if (myrank == root) {
    printf("Eigenvalue decomposition of Laplacian matrix\n");
    printf("solver = %s\n", solver_name);
    printf("dimension = %d\n", dim);
  }

  struct rokko_parallel_sparse_solver solver;
  rokko_parallel_sparse_solver_construct(&solver, solver_name, argc, argv);

  struct rokko_distributed_crs_matrix mat, Z;
  struct rokko_localized_vector w;

  rokko_distributed_crs_matrix_construct(&mat, dim, dim, solver);
  rokko_distributed_crs_matrix_construct(&Z, dim, dim, solver);

  int row;
  int row_start = rokko_distributed_crs_matrix_start_row(&mat);
  int row_end = rokko_distributed_crs_matrix_end_row(&mat);
  int cols[3];
  double values[3];
 
  if (row_start == 0) {
    values[0] = 1.;  values[1] = -1.;
    cols[0] = 0;   cols[1] = 1;
    ++row_start;
  }
  rokko_distributed_crs_matrix_insert(&mat, 0, 2, cols, values);

  if (row_end == (dim-1)) {
    --row_end;
  }
  
  values[0] = -1.;  values[1] = 2.;  values[2] = -1.;
  for (row = row_start; row < row_end; ++row) {
    cols[0] = row-1;   cols[1] = row;   cols[2] = row+1;
    rokko_distributed_crs_matrix_insert(&mat, row, 3, cols, values);
  }

  if (row_end == (dim-1)) {
    values[0] = -1.;  values[1] = 2.;
    cols[0] = dim-2;   cols[1] = dim-1;
  }
  rokko_distributed_crs_matrix_insert(&mat, dim-1, 2, cols, values);

  rokko_distributed_crs_matrix_complete(&mat);

  int nev = 10;
  int block_size = 5;
  int max_iters = 500;
  double tol = 1.0e-8;
  rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(&solver, &mat, nev, block_size, max_iters, tol);
  int num_conv = rokko_parallel_sparse_solver_num_conv(&solver);

  double eig_val = rokko_parallel_sparse_solver_eigenvalue(&solver, 0);
  int num_local_rows = rokko_distributed_crs_matrix_num_local_rows(&mat);
  double eig_vec[num_local_rows];
  int i = 0;
  rokko_parallel_sparse_solver_eigenvector(&solver, i, eig_vec);

  if (myrank == root) {
    printf("number of converged eigenpairs=%d\n", num_conv);
    printf("Computed Eigenvalue =\n");
    printf("%30.20f\n", eig_val);
    printf("Computed Eigenvector =\n");
    for (int j = 0; j < num_local_rows; ++j)
      printf("%30.20f ", eig_vec[j]);
    printf("\n");    
  }

  rokko_distributed_crs_matrix_destruct(&mat);
  rokko_distributed_crs_matrix_destruct(&Z);
  rokko_parallel_sparse_solver_destruct(&solver);

  MPI_Finalize();
  return 0;
}
