/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[]) {
  struct rokko_eigen_matrix mat;
  double *array_ptr;
  int dim, i, j;

  dim = 10;
  rokko_eigen_matrix_construct(&mat, dim, dim, rokko_matrix_col_major);
  array_ptr = rokko_eigen_matrix_get_array_pointer(mat);

  /* generate minij matrix */
  for(i=0; i<dim; ++i) {
    for(j=0; j<dim; ++j) {
      array_ptr[i + j*dim] = MIN(i, j) + 1;
    }
  }
  rokko_eigen_matrix_print(mat);

  rokko_eigen_matrix_destruct(&mat);

  return 0;
}
