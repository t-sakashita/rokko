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
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  int i, n;
  char** solver_names;

  n = rokko_serial_dense_ev_num_solvers();
  solver_names = rokko_serial_dense_ev_solvers();
  struct rokko_serial_dense_ev solver_sd;
  for (i = 0; i < n; ++i) {
    printf("%s\n", solver_names[i]);
    rokko_serial_dense_ev_construct(&solver_sd, solver_names[i], argc, argv);
    rokko_serial_dense_ev_destruct(&solver_sd);
  }
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);

  n = rokko_parallel_dense_ev_num_solvers();
  solver_names = rokko_parallel_dense_ev_solvers();
  struct rokko_parallel_dense_ev solver_pd;
  for (i = 0; i < n; ++i) {
    printf("%s\n", solver_names[i]);
    rokko_parallel_dense_ev_construct(&solver_pd, solver_names[i], argc, argv);
    rokko_parallel_dense_ev_destruct(&solver_pd);
  }
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);

  n = rokko_parallel_sparse_solver_num_solvers();
  solver_names = rokko_parallel_sparse_solver_solvers();
  struct rokko_parallel_sparse_solver solver_ps;
  for (i = 0; i < n; ++i) {
    printf("%s\n", solver_names[i]);
    rokko_parallel_sparse_solver_construct(&solver_ps, solver_names[i], argc, argv);
    rokko_parallel_sparse_solver_destruct(&solver_ps);
  }
  for (i = 0; i < n; ++i) free(solver_names[i]);
  free(solver_names);
}
