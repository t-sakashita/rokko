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
  printf("myrank=%d start_row=%d end_row=%d\n", vars->myrank, vars->start_row, vars->end_row);
  printf("myrank=%d num_local_rows=%d\n", vars->myrank, vars->num_local_rows);
}

void laplacian_multiply(const double* x, double* y, void* vars) {
  struct laplacian_vars* p = (struct laplacian_vars*)vars;
  int oo = p->dim;
  
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
  for (int k=1; k<p->end_k; ++k) {
    y[k] = - x[k-1] + 2 * x[k] - x[k+1];
  }
}

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
    printf("solver name = %s\n", solver_name);
    printf("Eigenvalue decomposition of Laplacian matrix\n");
    printf("solver = %s", solver_name);
    printf("dimension = %d\n", dim);
  }

  struct rokko_parallel_sparse_solver solver;
  rokko_parallel_sparse_solver_construct(&solver, solver_name, argc, argv);

  struct rokko_distributed_mfree mat;
  struct laplacian_vars vars;
  laplacian_initialize(dim, &vars);
  rokko_distributed_mfree_construct(&mat, laplacian_multiply, &vars, vars.dim, vars.num_local_rows);

  struct rokko_distributed_crs_matrix Z;
  rokko_distributed_crs_matrix_construct(&Z, dim, dim, solver);
  struct rokko_localized_vector w;

  int nev = 10;
  int block_size = 5;
  int max_iters = 500;
  double tol = 1.0e-8;
  rokko_parallel_sparse_solver_diagonalize_distributed_mfree(&solver, &mat, nev, block_size, max_iters, tol);
  int num_conv = rokko_parallel_sparse_solver_num_conv(&solver);

  double eig_val = rokko_parallel_sparse_solver_eigenvalue(&solver, 0);
  int num_local_rows = rokko_distributed_mfree_num_local_rows(&mat);
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

  rokko_distributed_mfree_destruct(&mat);
  rokko_distributed_crs_matrix_destruct(&Z);
  rokko_parallel_sparse_solver_destruct(&solver);

  MPI_Finalize();
  return 0;
}
