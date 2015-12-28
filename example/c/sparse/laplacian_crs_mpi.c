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

#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int num_solvers;
  char** solvers;
  if (argc >= 2) {
    num_solvers = 1;
    solvers = (char**)malloc((size_t)(num_solvers * sizeof(char*)));
    solvers[0] = (char*)malloc((size_t)((strlen(argv[1]) + 1) * sizeof(char)));
    strcpy(solvers[0], argv[1]);
  } else {
    num_solvers = rokko_parallel_sparse_ev_num_solvers();
    solvers = rokko_parallel_sparse_ev_solvers();
  }

  int dim = (argc == 3) ? dim = atoi(argv[2]) : 100;
  //  int nev = 10;
  //  int block_size = 5;
  //  int max_iters = 500;
  //  double tol = 1.0e-8;
  int s;
  for (s = 0; s < num_solvers; ++s) {
    struct rokko_parallel_sparse_ev solver;
    rokko_parallel_sparse_ev_construct(&solver, solvers[s], argc, argv);
    struct rokko_distributed_crs_matrix mat;
    rokko_distributed_crs_matrix_construct(&mat, dim, dim, solver);
    int row;
    int row_start = rokko_distributed_crs_matrix_start_row(&mat);
    int row_end = rokko_distributed_crs_matrix_end_row(&mat);
    int cols[3];
    double values[3];
    if (row_start == 0) {
      values[0] = 1.;  values[1] = -1.;
      cols[0] = 0;   cols[1] = 1;
      ++row_start;
      rokko_distributed_crs_matrix_insert(&mat, 0, 2, cols, values);
    }
    if (row_end == (dim-1)) {
      --row_end;
    }
    values[0] = -1.;  values[1] = 2.;  values[2] = -1.;
    for (row = row_start; row <= row_end; ++row) {
      cols[0] = row-1;   cols[1] = row;   cols[2] = row+1;
      rokko_distributed_crs_matrix_insert(&mat, row, 3, cols, values);
    }
    if (rokko_distributed_crs_matrix_end_row(&mat) == (dim-1)) {
      values[0] = -1.;  values[1] = 2.;
      cols[0] = dim-2;   cols[1] = dim-1;
      rokko_distributed_crs_matrix_insert(&mat, dim-1, 2, cols, values);
    }
    rokko_distributed_crs_matrix_complete(&mat);
    if (rank == 0) {
      printf("Eigenvalue decomposition of Laplacian matrix\n");
      printf("solver = %s\n", solvers[s]);
      printf("dimension = %d\n", dim);
    }

    struct rokko_parameters params;
    rokko_parameters_construct(&params);
    // set some parameters
    rokko_parameters_set_int(&params, "max_block_size", 5);
    rokko_parameters_set_int(&params, "max_iters", 500);
    rokko_parameters_set_double(&params, "conv_tol", 1.0e-12);
    rokko_parameters_set_int(&params, "num_eigvals", 1);
    rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(solver, mat, params);

    int num_conv = rokko_parallel_sparse_ev_num_conv(&solver);
    if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
    int num_local_rows = rokko_distributed_crs_matrix_num_local_rows(&mat);
    double eig_vec[num_local_rows];
    int i, j;
    if (rank == 0) {
      printf("number of converged eigenpairs = %d\n", num_conv);
      printf("smallest eigenvalues: ");
      for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_parallel_sparse_ev_eigenvalue(&solver, i));
      printf("\n");
      rokko_parallel_sparse_ev_eigenvector(&solver, 0, eig_vec);
      printf("smallest eigenvector: ");
      for (j = 0; j < num_local_rows; ++j)
        printf("%30.20f ", eig_vec[j]);
      printf("\n");    
    }

    rokko_distributed_crs_matrix_destruct(&mat);
    rokko_parallel_sparse_ev_destruct(&solver);
  }
  
  for (s = 0; s < num_solvers; ++s) free(solvers[s]);
  free(solvers);
  MPI_Finalize();
  return 0;
}
