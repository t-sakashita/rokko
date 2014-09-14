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
#include <rokko/rokko_sparse.h>

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int provided, ierr, myrank, nprocs;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int root = 0;

  int nev = 1;
  int blockSize = 1;
  int maxIters = 500;
  double tol = 1.0e-8;

  int L = 3;
  int dim = 1 << L;
  int lattice_first[L], lattice_second[L];
  int count, l;
  count = 0;
  for (l = 0; l < L; ++l) {
    lattice_first[count] = l;
    lattice_second[count] = (l+1) % L;
  }

  char solver_name[7] = "anasazi";

  struct rokko_distributed_crs_matrix mat, Z;
  struct rokko_parallel_sparse_solver solver;

  struct rokko_localized_vector w;

  printf("solver name = %s\n", solver_name);
  printf("matrix dimension = %d\n", dim);

  rokko_parallel_sparse_solver_construct(&solver, solver_name, argc, argv);
  rokko_distributed_crs_matrix_construct(&mat, dim, dim, solver);
  rokko_distributed_crs_matrix_construct(&Z, dim, dim, solver);
  rokko_localized_vector_construct(&w, dim);

  int row;
  int row_start = rokko_distributed_crs_matrix_start_row(&mat);
  int row_end = rokko_distributed_crs_matrix_end_row(&mat);
  int cols[dim];
  double values[dim];
 
  int diag;
  int i, j, m1, m2, m3;
  count = 0;
  for (row = row_start; row < row_end; ++row) {
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
    rokko_distributed_crs_matrix_insert(&mat, count, row, cols, values);
  }

  rokko_distributed_crs_matrix_complete(&mat);
  //rokko_solver_diagonalize_distributed_crs_matrix(&solver, &mat);

  /*
  rokko_distributed_crs_matrix_print(mat);

  if (myrank == 0) {
    printf("Computed Eigenvalues =\n");
    for (i = 0; i < dim; ++i)
      printf("%30.20f\n", rokko_localized_vector_get(w, i));
  }
  */
  rokko_distributed_crs_matrix_destruct(&mat);
  rokko_distributed_crs_matrix_destruct(&Z);
  rokko_localized_vector_destruct(&w);

  rokko_parallel_sparse_solver_destruct(&solver);

  MPI_Finalize();
  return 0;
}
