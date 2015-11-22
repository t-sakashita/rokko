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

  int L = (argc >= 3) ? atoi(argv[2]) : 10;
  int dim = 1 << L;
  int lattice_first[L], lattice_second[L];
  int l;
  for (l = 0; l < L; ++l) {
    lattice_first[l] = l;
    lattice_second[l] = (l+1) % L;
  }

  int nev = 10;
  int block_size = 5;
  int max_iters = 500;
  double tol = 1.0e-8;
  int s;
  for (s = 0; s < num_solvers; ++s) {
    struct rokko_parallel_sparse_ev solver;
    rokko_parallel_sparse_ev_construct(&solver, solvers[s], argc, argv);
    struct rokko_distributed_crs_matrix mat;
    rokko_distributed_crs_matrix_construct(&mat, dim, dim, solver);
    int row;
    int row_start = rokko_distributed_crs_matrix_start_row(&mat);
    int row_end = rokko_distributed_crs_matrix_end_row(&mat);
    int cols[dim];
    double values[dim];
    int count;
    double diag;
    int i, j, m1, m2, m3;
    for (row = row_start; row <= row_end; ++row) {
      count = 0;
      diag = 0;
      for (l = 0; l < L; ++l) {
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
    if (rank == 0) {
      printf("Eigenvalue decomposition of antiferromagnetic Heisenberg chain\n");
      printf("solver = %s\n", solvers[s]);
      printf("L = %d\n", L);
      printf("dimension = %d\n", dim);
    }

    rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(&solver, &mat, nev, block_size, max_iters, tol);

    int num_conv = rokko_parallel_sparse_ev_num_conv(&solver);
    if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
    int num_local_rows = rokko_distributed_crs_matrix_num_local_rows(&mat);
    double eig_vec[num_local_rows];
    rokko_parallel_sparse_ev_eigenvector(&solver, 0, eig_vec);
    if (rank == 0) {
      printf("number of converged eigenpairs = %d\n", num_conv);
      printf("smallest eigenvalues: ");
      for (i = 0; i < num_conv; ++i) printf("%30.20f", rokko_parallel_sparse_ev_eigenvalue(&solver, i));
      printf("\n");
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
