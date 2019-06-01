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

#include <cblas.h>
#include <lapacke.h>
#include <cmatrix.h>
#include <math.h>

int imin(int x, int y) { return (x < y) ? x : y; }

int main(int argc, char** argv) {
  int n = 5;
  int i, j, info;
  float norm, check;
  float **a;
  float *work;

  if (argc > 1) n = atoi(argv[1]);

  // generate matrix and rhs vector
  a = alloc_smatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      mat_elem(a, i, j) = imin(i+1, j+1);
    }
  }
  printf("Matrix A: ");
  fprint_smatrix(stdout, n, n, a);
  work = alloc_svector(n);

  /* calculate various norms */
  norm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'M', n, n, mat_ptr(a), n, vec_ptr(work));
  printf("element of largest absolute value = %e\n", norm);

  norm = LAPACKE_slange(LAPACK_COL_MAJOR, '1', n, n, mat_ptr(a), n);
  printf("one norm = %e\n", norm);

  norm = LAPACKE_slange(LAPACK_COL_MAJOR, 'I', n, n, mat_ptr(a), n);
  printf("infinity norm = %e\n", norm);

  norm = LAPACKE_slange(LAPACK_COL_MAJOR, 'F', n, n, mat_ptr(a), n);
  printf("Frobenius norm = %e\n", norm);

  /* check of result */
  check = 0;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      check += mat_elem(a, i, j) * mat_elem(a, i, j);
    }
  }
  if (fabsf(check - norm * norm) > 1e-6 * check) {
    fprintf(stderr, "Error: check error %e != %e\n", norm * norm, check);
    exit(255);
  }

  free_smatrix(a);
  return 0;
}
