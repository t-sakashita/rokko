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

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define DIM 10

int main(int argc, char *argv[]) {
  struct rokko_eigen_matrix mat;
  double array_ptr[DIM*DIM];
  int i, j;

  rokko_eigen_matrix_construct_array_sizes(&mat, DIM, DIM, array_ptr, rokko_matrix_col_major);

  /* generate frank matrix */
  for(i=0; i<DIM; ++i) {
    for(j=0; j<DIM; ++j) {
      array_ptr[i + j*DIM] = DIM - MAX(i, j);
    }
  }
  rokko_eigen_matrix_print(mat);

  /*rokko_eigen_matrix_destruct(&mat);*/

  return 0;
}
