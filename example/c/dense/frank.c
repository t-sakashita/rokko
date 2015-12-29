/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*    
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.h>
#include <rokko/utility/frank_matrix_c.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int dim;
  struct rokko_serial_dense_ev solver;
  struct rokko_localized_matrix mat, Z;
  struct rokko_localized_vector w;
  char* solver_name;

  if (argc == 3) {
    solver_name = argv[1];
    dim = atoi(argv[2]);
  } else {
    fprintf(stderr, "error: %s solver_name dimension\n", argv[0]);
    exit(127);
  }
    
  printf("solver name = %s\n", solver_name);
  printf("matrix dimension = %d\n", dim);

  rokko_serial_dense_ev_construct(&solver, solver_name, argc, argv);

  rokko_localized_matrix_construct(&mat, dim, dim, rokko_matrix_col_major);
  rokko_localized_matrix_construct(&Z, dim, dim, rokko_matrix_col_major);
  rokko_localized_vector_construct(&w, dim);

  /* generate frank matrix */
  rokko_frank_matrix_generate_localized_matrix(mat);
  rokko_localized_matrix_print(mat);

  rokko_serial_dense_ev_diagonalize_localized_matrix(solver, mat, w, Z);

  printf("Computed Eigenvalues =\n");
  int i;
  for (i = 0; i < dim; ++i)
    printf("%30.20f\n", rokko_localized_vector_get(w, i));

  rokko_localized_matrix_destruct(&mat);
  rokko_localized_matrix_destruct(&Z);
  rokko_localized_vector_destruct(&w);
  rokko_serial_dense_ev_destruct(&solver);

  return 0;
}
