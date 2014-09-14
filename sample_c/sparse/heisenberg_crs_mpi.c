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
  int provided, ierr, myrank, nprocs;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int root = 0;

  int L = 8;
  int dim = 1 << L;
  int lattice_first[L], lattice_second[L];
  int l;
  for (l = 0; l < L; ++l) {
    lattice_first[l] = l;
    lattice_second[l] = (l+1) % L;
  }

  char solver_name[7] = "anasazi";
  printf("solver name = %s\n", solver_name);
  if (myrank == root) {
    printf("Eigenvalue decomposition of antiferromagnetic Heisenberg chain\n");
    printf("solver = anasazi/Krylov-Schur\n");
    printf("L = %d\n", L);
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
  int cols[dim];
  double values[dim];
 
  int count;
  double diag;
  int i, j, m1, m2, m3;
  for (row = row_start; row < row_end; ++row) {
    count = 0;
    diag = 0;
    for (l = 0;  l < L; ++l) {
      i = lattice_first[l];
      j = lattice_second[l];
      m1 = 1 << i;
      m2 = 1 << j;
      m3 = m1 + m2;
      if (((row & m3) == m1) || ((row & m3) == m2)) {
        cols[count] = row^m3;
        values[count] = 0.5;
	++count;
        diag += -0.25;
      } else {
        diag += 0.25;
      }
    }
    cols[count] = row;
    values[count] = diag;
    ++count;
    rokko_distributed_crs_matrix_insert(&mat, row, count, cols, values);
  }
  rokko_distributed_crs_matrix_complete(&mat);

  int nev = 10;
  int blockSize = 5;
  int maxIters = 500;
  double tol = 1.0e-8;
  rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(&solver, &mat, nev, blockSize, maxIters, tol);
  int num_conv = rokko_parallel_sparse_solver_num_conv(&solver);

  double eig_val = rokko_parallel_sparse_solver_eigenvalue(&solver, 0);
  //double* eig_vec;
  //eig_vec = (double *)malloc(sizeof(double) * num_local_rows);
  int num_local_rows = rokko_distributed_crs_matrix_num_local_rows(&mat);
  double eig_vec[num_local_rows];
  rokko_parallel_sparse_solver_eigenvector(&solver, i, eig_vec);

  if (myrank == root) {
    printf("number of converged eigenpairs=%d\n", num_conv);
    printf("Computed Eigenvalue =\n");
    printf("%30.20f\n", eig_val);
    printf("Computed Eigenvector =\n");
    for (j = 0; j < num_local_rows; ++j)
      printf("%30.20f ", eig_vec[j]);
    printf("\n");    
  }

  rokko_distributed_crs_matrix_destruct(&mat);
  rokko_distributed_crs_matrix_destruct(&Z);
  rokko_parallel_sparse_solver_destruct(&solver);

  MPI_Finalize();
  return 0;
}
