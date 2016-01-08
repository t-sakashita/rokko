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
  struct rokko_parameters params;
  char* library_routine, *library, *routine;

  library_routine = rokko_serial_dense_ev_default_solver();
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = atoi(argv[2]);
  rokko_split_solver_name(library_routine, &library, &routine);
  rokko_serial_dense_ev_construct(&solver, library, argc, argv);

  printf("library = %s\n", library);
  printf("routine = %s\n", routine);
  printf("dimension = %d\n", dim);
  
  rokko_localized_matrix_construct(&mat, dim, dim, rokko_matrix_col_major);
  rokko_localized_matrix_construct(&Z, dim, dim, rokko_matrix_col_major);
  rokko_localized_vector_construct(&w, dim);
  rokko_parameters_construct(&params);
  rokko_parameters_set_string(params, "routine", routine);
  
  /* generate frank matrix */
  rokko_frank_matrix_generate_localized_matrix(mat);
  rokko_localized_matrix_print(mat);

  rokko_serial_dense_ev_diagonalize_params(solver, mat, w, Z, params);

  printf("Computed Eigenvalues =\n");
  int i;
  for (i = 0; i < dim; ++i)
    printf("%30.20f\n", rokko_localized_vector_get(w, i));

  rokko_parameters_destruct(&params);
  rokko_localized_matrix_destruct(&mat);
  rokko_localized_matrix_destruct(&Z);
  rokko_localized_vector_destruct(&w);
  rokko_serial_dense_ev_destruct(&solver);

  return 0;
}
