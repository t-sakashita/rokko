/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.h>
#include <stdio.h>

int main(/* int argc , char** argv */) {
  printf("%s\n", rokko_serial_dense_ev_default_solver());

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  printf("%s\n", rokko_parallel_dense_ev_default_solver());
#endif

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  printf("%s\n", rokko_parallel_sparse_ev_default_solver());
#endif

  return 0;
}
