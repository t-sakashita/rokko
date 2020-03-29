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
#include <stdlib.h>

int main(/* int argc , char** argv */) {
  int i, n;
  char** solver_names;

  n = rokko_serial_dense_ev_num_solvers();
  solver_names = rokko_serial_dense_ev_solvers();
  for (i = 0; i < n; ++i)
    printf("%s\n", solver_names[i]);
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  n = rokko_parallel_dense_ev_num_solvers();
  solver_names = rokko_parallel_dense_ev_solvers();
  for (i = 0; i < n; ++i)
    printf("%s\n", solver_names[i]);
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);
#endif

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  n = rokko_parallel_sparse_ev_num_solvers();
  solver_names = rokko_parallel_sparse_ev_solvers();
  for (i = 0; i < n; ++i)
    printf("%s\n", solver_names[i]);
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);
#endif

  return 0;
}
