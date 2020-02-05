/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
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

#define min(a, b) ((a) < (b) ? (a) : (b))

struct laplacian_vars {
  MPI_Comm comm;
  int nprocs, myrank;
  int dim, local_offset, num_local_rows;
  int start_row, end_row;
  int end_k;
  int is_first_proc, is_last_proc;
} lap_var;


void laplacian_multiply(const double *const x, double *const y) {
  double buf_m, buf_p;
  MPI_Status status_m, status_p;

  if (lap_var.num_local_rows == 0) return;

  if ((!lap_var.is_first_proc) && (lap_var.nprocs != 1)) {
    //std::cout << "recv myrank=" << lap_var.myrank << std::endl;
    MPI_Send(&x[0], 1, MPI_DOUBLE, lap_var.myrank-1, 0, lap_var.comm);
    MPI_Recv(&buf_m, 1, MPI_DOUBLE, lap_var.myrank-1, 0, lap_var.comm, &status_m);
    //std::cout << "buffff=" << buf << std::endl;
  }

  if ((!lap_var.is_last_proc) && (lap_var.nprocs != 1)) {
    //std::cout << "send myrank=" << lap_var.myrank << std::endl;
    MPI_Recv(&buf_p, 1, MPI_DOUBLE, lap_var.myrank+1, 0, lap_var.comm, &status_p);
    MPI_Send(&x[lap_var.end_k], 1, MPI_DOUBLE, lap_var.myrank+1, 0, lap_var.comm);
    //std::cout << "buffff=" << buf2 << std::endl;
  }

  if (lap_var.is_first_proc) {
    if (lap_var.num_local_rows != 1) {
      y[0] = x[0] - x[1];
      if (lap_var.nprocs != 1) y[lap_var.end_k] = - x[lap_var.end_k - 1] + 2 * x[lap_var.end_k] - buf_p;
    }
    else {
      y[0] = x[0] - buf_p;
    }
  }

  if (lap_var.is_last_proc) {
    if (lap_var.num_local_rows != 1) {
      if (lap_var.nprocs != 1) y[0] = - buf_m + 2 * x[0] - x[1];
      y[lap_var.end_k] = 2 * x[lap_var.end_k] - x[lap_var.end_k - 1];
      }
    else {
      y[lap_var.end_k] = 2 * x[lap_var.end_k] - buf_m;
    }
  }
  if (!(lap_var.is_first_proc || lap_var.is_last_proc)) { // neither first or last process
    if (lap_var.num_local_rows != 1) {
      y[0] = - buf_m + 2 * x[0] - x[1];
      y[lap_var.end_k] = - x[lap_var.end_k - 1] + 2 * x[lap_var.end_k] - buf_p;
    }
    else {
      y[0] = - buf_m + 2 * x[0] - buf_p;
    }
  }
  // from 1 to end-1
  int k;
  for (k=1; k<lap_var.end_k; ++k) {
    y[k] = - x[k-1] + 2 * x[k] - x[k+1];
  }
}

void laplacian_initialize(struct rokko_distributed_mfree* mat, int dim, struct laplacian_vars* vars) {
  rokko_distributed_mfree_construct_no_context(mat, laplacian_multiply, dim, MPI_COMM_WORLD);

  lap_var.comm = MPI_COMM_WORLD;
  MPI_Comm_size(lap_var.comm, &lap_var.nprocs);
  MPI_Comm_rank(lap_var.comm, &lap_var.myrank);

  lap_var.dim = dim;
  lap_var.num_local_rows = rokko_distributed_mfree_num_local_rows(*mat);
  lap_var.start_row = rokko_distributed_mfree_start_row(*mat);
  lap_var.end_row = rokko_distributed_mfree_end_row(*mat);

  lap_var.is_first_proc = (lap_var.start_row == 0);
  lap_var.is_last_proc = (lap_var.end_row == (lap_var.dim));

  lap_var.end_k = lap_var.num_local_rows - 1;
}

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
  struct rokko_distributed_mfree mat;
  struct laplacian_vars vars;
  laplacian_initialize(&mat, dim, &vars);
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
  rokko_parallel_sparse_ev_diagonalize_distributed_mfree(solver, mat, params);

  int num_conv = rokko_parallel_sparse_ev_num_conv(solver);
  if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
  int num_local_rows = rokko_distributed_mfree_num_local_rows(mat);
  double eig_vec[num_local_rows];
  rokko_parallel_sparse_ev_eigenvector(solver, 0, eig_vec);
  if (rank == 0) {
    printf("number of converged eigenpairs = %d\n", num_conv);
    printf("largest eigenvalues: ");
    int i, j;
    for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_parallel_sparse_ev_eigenvalue(solver, i));
    printf("\n");
    printf("smallest theoretical eigenvalues: ");
    for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_laplacian_matrix_eigenvalue(dim, i));
    //for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_laplacian_matrix_eigenvalue(dim, dim-1-i));
    printf("\n");
    printf("largest eigenvector: ");
    for (j = 0; j < num_local_rows; ++j)
      printf("%30.20f ", eig_vec[j]);
    printf("\n");
  }
  rokko_distributed_mfree_destruct(&mat);
  rokko_parallel_sparse_ev_destruct(&solver);

  MPI_Finalize();
  return 0;
}
