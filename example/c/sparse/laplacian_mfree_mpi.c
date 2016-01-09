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

#define min(a, b) ((a) < (b) ? (a) : (b))

struct laplacian_vars {
  MPI_Comm comm;
  int nprocs, myrank;
  int dim, local_offset, num_local_rows;
  int start_row, end_row;
  int start_k, end_k;
  int is_first_proc, is_last_proc;
  double buf_m, buf_p;
  MPI_Status status_m, status_p;
};

void laplacian_initialize(int dim, struct laplacian_vars* vars) {
  vars->comm = MPI_COMM_WORLD;
  MPI_Comm_size(vars->comm, &vars->nprocs);
  MPI_Comm_rank(vars->comm, &vars->myrank);

  vars->dim = dim;
  int tmp = vars->dim / vars->nprocs;
  int rem = vars->dim % vars->nprocs;
  vars->num_local_rows = (vars->dim + vars->nprocs - vars->myrank - 1) / vars->nprocs;
  vars->start_row = tmp * vars->myrank + min(rem, vars->myrank);
  vars->end_row = vars->start_row + vars->num_local_rows - 1;

  if (vars->start_row == 0)  vars->is_first_proc = 1;
  else vars->is_first_proc = 0;
  
  if (vars->end_row == (vars->dim-1))  vars->is_last_proc = 1;
  else vars->is_last_proc = 0;
    
  vars->end_k = vars->num_local_rows - 1;
}

void laplacian_multiply(const double* x, double* y, void* vars) {
  struct laplacian_vars* p = (struct laplacian_vars*)vars;
  
  if (p->num_local_rows == 0) return;
    
  if ((!p->is_first_proc) && (p->nprocs != 1)) {
    //std::cout << "recv myrank=" << p->myrank << std::endl;
    MPI_Send(&x[0], 1, MPI_DOUBLE, p->myrank-1, 0, p->comm);
    MPI_Recv(&p->buf_m, 1, MPI_DOUBLE, p->myrank-1, 0, p->comm, &p->status_m);
    //std::cout << "buffff=" << buf << std::endl;
  }
  
  if ((!p->is_last_proc) && (p->nprocs != 1)) {
    //std::cout << "send myrank=" << p->myrank << std::endl;
    MPI_Recv(&p->buf_p, 1, MPI_DOUBLE, p->myrank+1, 0, p->comm, &p->status_p);
    MPI_Send(&x[p->end_k], 1, MPI_DOUBLE, p->myrank+1, 0, p->comm);
    //std::cout << "buffff=" << buf2 << std::endl;
  }

  if (p->is_first_proc) {
    if (p->num_local_rows != 1) {
      y[0] = x[0] - x[1];
      if (p->nprocs != 1) y[p->end_k] = - x[p->end_k - 1] + 2 * x[p->end_k] - p->buf_p;
    }
    else {
      y[0] = x[0] - p->buf_p;
    }
  }
  
  if (p->is_last_proc) {
    if (p->num_local_rows != 1) {
      if (p->nprocs != 1) y[0] = - p->buf_m + 2 * x[0] - x[1];
      y[p->end_k] = 2 * x[p->end_k] - x[p->end_k - 1];
      }
    else {
      y[p->end_k] = 2 * x[p->end_k] - p->buf_m;
    }
  }
  if (!(p->is_first_proc || p->is_last_proc)) { // neither first or last process
    if (p->num_local_rows != 1) {
      y[0] = - p->buf_m + 2 * x[0] - x[1];
      y[p->end_k] = - x[p->end_k - 1] + 2 * x[p->end_k] - p->buf_p;
    }
    else {
      y[0] = - p->buf_m + 2 * x[0] - p->buf_p;
    }
  }
  // from 1 to end-1
  int k;
  for (k=1; k<p->end_k; ++k) {
    y[k] = - x[k-1] + 2 * x[k] - x[k+1];
  }
}

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

  int s;
  for (s = 0; s < num_solvers; ++s) {
    struct rokko_parallel_sparse_ev solver;
    rokko_parallel_sparse_ev_construct(&solver, solvers[s], argc, argv);
    struct rokko_distributed_mfree mat;
    struct laplacian_vars vars;
    laplacian_initialize(dim, &vars);
    rokko_distributed_mfree_construct(&mat, laplacian_multiply, &vars, vars.dim, vars.num_local_rows);
    if (rank == 0) {
      printf("Eigenvalue decomposition of Laplacian matrix\n");
      printf("solver = %s\n", solvers[s]);
      printf("dimension = %d\n", dim);
    }

    struct rokko_parameters params;
    rokko_parameters_construct(&params);
    // set some parameters
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
      printf("largest eigenvector: ");
      for (j = 0; j < num_local_rows; ++j)
        printf("%30.20f ", eig_vec[j]);
      printf("\n");    
    }
    rokko_distributed_mfree_destruct(&mat);
    rokko_parallel_sparse_ev_destruct(&solver);
  }

  for (s = 0; s < num_solvers; ++s) free(solvers[s]);
  free(solvers);
  MPI_Finalize();
  return 0;
}
