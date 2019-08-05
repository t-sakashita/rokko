/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <rokko/cmatrix.h>

int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int n = 5;
  int i, j, info;
  double norm, check;
  double **a;

  if (argc > 1) n = atoi(argv[1]);

  // generate matrix and rhs vector
  a = alloc_dmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      mat_elem(a, i, j) = n - 0.253 * imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix(stdout, n, n, a);

  /* calculate various norms */
  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, mat_ptr(a), n);
  printf("norm1(A) = %e\n", norm);

  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', n, n, mat_ptr(a), n);
  printf("normI(A) = %e\n", norm);

  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', n, n, mat_ptr(a), n);
  printf("normF(A) = %e\n", norm);

  /* check of result */
  check = 0;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      check += mat_elem(a, i, j) * mat_elem(a, i, j);
    }
  }
  if (fabs(check - norm * norm) > 1e-12 * check) {
    fprintf(stderr, "Error: check error %e != %e\n", norm * norm, check);
    exit(255);
  }

  free_dmatrix(a);
  return 0;
}
