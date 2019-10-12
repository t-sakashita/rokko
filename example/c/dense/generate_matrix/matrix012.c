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
#include <rokko/utility/matrix012.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  struct rokko_eigen_matrix mat;
  int dim;

  dim = 10;
  rokko_eigen_matrix_construct(&mat, dim, dim, rokko_matrix_col_major);
  rokko_matrix012_generate_eigen_matrix(mat);
  rokko_eigen_matrix_print(mat);

  rokko_eigen_matrix_destruct(&mat);

  return 0;
}
