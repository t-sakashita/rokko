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
#include <rokko/utility/frank_matrix.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int dim = 10;
  struct rokko_parallel_dense_ev solver;
  struct rokko_grid grid;
  struct rokko_mapping_bc map;
  struct rokko_distributed_matrix mat, Z;
  struct rokko_eigen_vector w;
  struct rokko_parameters params;
  char *library_routine, *library, *routine;

  int provided, myrank, nprocs, i;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (argc >= 2) library_routine = argv[1];
  else library_routine = rokko_parallel_dense_ev_default_solver();
  if (argc >= 3) dim = atoi(argv[2]);
  rokko_split_solver_name(library_routine, &library, &routine);

  if (myrank == 0) {
    printf("library = %s\n", library);
    printf("routine = %s\n", routine);
    printf("dimension = %d\n", dim);
  }

  rokko_parallel_dense_ev_construct(&solver, library, argc, argv);
  rokko_grid_construct(&grid, MPI_COMM_WORLD, rokko_grid_row_major);
  map = rokko_parallel_dense_ev_default_mapping(solver, dim, grid);
  rokko_distributed_matrix_construct(&mat, map);
  rokko_distributed_matrix_construct(&Z, map);
  rokko_eigen_vector_construct(&w, dim);

  /* generate frank matrix */
  rokko_frank_matrix_generate_distributed_matrix(mat);
  rokko_distributed_matrix_print(mat);

  rokko_parameters_construct(&params);
  rokko_parameters_set_string(params, "routine", routine);
  rokko_parallel_dense_ev_diagonalize(solver, mat, w, Z, params);

  if (myrank == 0) {
    printf("Computed Eigenvalues =\n");
    for (i = 0; i < dim; ++i)
      printf("%30.20f\n", rokko_eigen_vector_get(w, i));
  }

  rokko_distributed_matrix_destruct(&mat);
  rokko_distributed_matrix_destruct(&Z);
  rokko_eigen_vector_destruct(&w);
  rokko_parallel_dense_ev_destruct(&solver);
  rokko_mapping_bc_destruct(&map);
  rokko_grid_destruct(&grid);

  MPI_Finalize();
  return 0;
}
