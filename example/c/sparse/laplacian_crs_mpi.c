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

#include <rokko/utility/laplacian_matrix.h>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char* library_routine, *library, *routine;
  if (argc >= 2) {
    library_routine = rokko_parallel_sparse_ev_default_solver();
    if (argc >= 2) library_routine = argv[1];
  } else {
    library_routine = rokko_parallel_sparse_ev_default_solver();
  }
  rokko_split_solver_name(library_routine, &library, &routine);
  int dim = (argc == 3) ? dim = atoi(argv[2]) : 100;

  struct rokko_parallel_sparse_ev solver;
  rokko_parallel_sparse_ev_construct(&solver, library, argc, argv);
  struct rokko_distributed_crs_matrix mat;
  rokko_distributed_crs_matrix_construct(&mat, dim, dim, solver);
  int row;
  int row_start = rokko_distributed_crs_matrix_start_row(mat);
  int row_end = rokko_distributed_crs_matrix_end_row(mat);
  int cols[3];
  double values[3];
  printf("row_start=%d, row_end=%d\n", row_start, row_end);
  if (row_start == 0) {
    values[0] = 1.;  values[1] = -1.;
    cols[0] = 0;   cols[1] = 1;
    ++row_start;
    rokko_distributed_crs_matrix_insert(mat, 0, 2, cols, values);
  }
  int row_loop_end = row_end;
  if (row_end == dim) {
    --row_loop_end;
  }
  values[0] = -1.;  values[1] = 2.;  values[2] = -1.;
  for (row = row_start; row < row_loop_end; ++row) {
    cols[0] = row-1;   cols[1] = row;   cols[2] = row+1;
    rokko_distributed_crs_matrix_insert(mat, row, 3, cols, values);
  }
  if (row_end == dim) {
    values[0] = -1.;  values[1] = 2.;
    cols[0] = dim-2;   cols[1] = dim-1;
    rokko_distributed_crs_matrix_insert(mat, dim-1, 2, cols, values);
  }
  rokko_distributed_crs_matrix_complete(mat);
  if (rank == 0) {
    printf("Eigenvalue decomposition of Laplacian matrix\n");
    printf("solver = %s\n", library);
    printf("dimension = %d\n", dim);
  }
  
  struct rokko_parameters params;
  rokko_parameters_construct(&params);
  // set some parameters
  if (routine[0] != '\0')  rokko_parameters_set_string(params, "routine", routine);
  rokko_parameters_set_int(params, "block_size", 5);
  rokko_parameters_set_int(params, "max_iters", 500);
  rokko_parameters_set_double(params, "conv_tol", 1.0e-8);
  rokko_parameters_set_int(params, "num_eigvals", 1);
  rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(solver, mat, params);
  
  int num_conv = rokko_parallel_sparse_ev_num_conv(solver);
  if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
  int num_local_rows = rokko_distributed_crs_matrix_num_local_rows(mat);
  double eig_vec[num_local_rows];
  int i, j;
  if (rank == 0) {
    printf("number of converged eigenpairs = %d\n", num_conv);
    printf("smallest eigenvalues: ");
    for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_parallel_sparse_ev_eigenvalue(solver, i));
    printf("\n");
    printf("smallest theoretical eigenvalues: ");
    for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_laplacian_matrix_eigenvalue(dim, i));
    //for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_laplacian_matrix_eigenvalue(dim, dim-1-i));
    printf("\n");
    rokko_parallel_sparse_ev_eigenvector(solver, 0, eig_vec);
    printf("smallest eigenvector: ");
    for (j = 0; j < num_local_rows; ++j)
      printf("%30.20f ", eig_vec[j]);
    printf("\n");    
  }
  
  rokko_distributed_crs_matrix_destruct(&mat);
  rokko_parallel_sparse_ev_destruct(&solver);
  
  MPI_Finalize();
  return 0;
}
