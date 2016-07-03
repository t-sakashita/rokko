/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <rokko/rokko.h>
#include <rokko/utility/frank_matrix_c.h>
#include <stdio.h>
#include <stdlib.h>

void diagonalize_fixedB(rokko_parallel_dense_ev* solver, rokko::distributed_matrix<T, MATRIX_MAJOR>* A, rokko::distributed_matrix<T, MATRIX_MAJOR>& B, rokko::localized_vector<double>& eigval, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvec, T tol = 0) {
  rokko::distributed_matrix<double, matrix_major> tmp(A.get_mapping()), Binvroot(A.get_mapping()), mat(A.get_mapping());
  rokko::parameters params;
  int myrank = A.get_myrank();
  params.set("routine", "");
  solver.diagonalize(B, eigval, eigvec, params);
  // computation of B^{-1/2}
  for(int i=0; i<eigval.size(); ++i)
    eigval(i) = (eigval(i) > tol) ? sqrt(1/eigval(i)) : 0;
  function_matrix(eigval, eigvec, Binvroot, tmp);
  
  // computation of B^{-1/2} A B^{-1/2}
  product(1, Binvroot, false, A, false, 0, tmp);
  product(1, tmp, false, Binvroot, false, 0, mat);
  // diagonalization of B^{-1/2} A B^{-1/2}
  solver.diagonalize(mat, eigval, tmp, params);

  // computation of {eigvec of Ax=lambda Bx} = B^{-1/2} {eigvec of B^{-1/2} A B^{-1/2}}
  product(1, Binvroot, false, tmp, false, 0, eigvec);
}

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_parallel_dense_ev solver;
  struct rokko_distributed_matrix mat, Z;
  struct rokko_grid grid;
  struct rokko_localized_vector w;
  char* solver_name;

  int provided, ierr, myrank, nprocs, i;

  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (argc == 3) {
    solver_name = argv[1];
    dim = atoi(argv[2]);
  } else {
    fprintf(stderr, "error: %s solver_name dimension\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 34);
  }
    
  printf("solver name = %s\n", solver_name);
  printf("matrix dimension = %d\n", dim);

  rokko_parallel_dense_ev_construct(&solver, solver_name, argc, argv);
  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);

  rokko_distributed_matrix_construct(&mat, dim, dim, grid, solver, rokko_matrix_col_major);
  rokko_distributed_matrix_construct(&Z, dim, dim, grid, solver, rokko_matrix_col_major);
  rokko_localized_vector_construct(&w, dim);

  /* generate frank matrix */
  rokko_frank_matrix_generate_distributed_matrix(&mat);
  rokko_distributed_matrix_print(mat);

  rokko_parallel_dense_ev_diagonalize_distributed_matrix(&solver, &mat, &w, &Z);

  if (myrank == 0) {
    printf("Computed Eigenvalues =\n");
    for (i = 0; i < dim; ++i)
      printf("%30.20f\n", rokko_localized_vector_get(w, i));
  }

  rokko_distributed_matrix_destruct(&mat);
  rokko_distributed_matrix_destruct(&Z);
  rokko_localized_vector_destruct(&w);
  rokko_parallel_dense_ev_destruct(&solver);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
