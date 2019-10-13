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

#define DIM 5

int main(int argc, char *argv[]) {
  struct rokko_eigen_vector vector;
  double array[DIM];
  int i;

  rokko_eigen_vector_construct_by_array(&vector, DIM, array);

  for(i=0; i<DIM; ++i)
    array[i] = i + 1;

  rokko_eigen_vector_print(vector);

  /*rokko_eigen_vector_destruct(&vector);*/

  return 0;
}
