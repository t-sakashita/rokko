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

int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int n = 5;
  int i, j, info;
  double norm;
  double *x, *b;
  double **a, **lu;
  int *ipiv;

  if (argc > 1) n = atoi(argv[1]);

  // generate matrix and rhs vector
  a = alloc_dmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      mat_elem(a, i, j) = n - imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix(stdout, n, n, a);
  b = alloc_dvector(n);
  for (i = 0; i < n; ++i) b[i] = i * i + 1;
  printf("Vector b: ");
  fprint_dvector(stdout, n, b);

  /* solve linear equation */
  lu = alloc_dmatrix(n, n);
  cblas_dcopy(n * n, mat_ptr(a), 1, mat_ptr(lu), 1);
  x = alloc_dvector(n);
  cblas_dcopy(n, vec_ptr(b), 1, vec_ptr(x), 1);
  ipiv = alloc_ivector(n);
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, mat_ptr(lu), n, vec_ptr(ipiv));
  if (info != 0) {
    fprintf(stderr, "Error: dgetrf fails\n");
    exit(255);
  }
  info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'n', n, 1, mat_ptr(lu), n, vec_ptr(ipiv),
                        vec_ptr(x), n);
  if (info != 0) {
    fprintf(stderr, "Error: dgetrs fails\n");
    exit(255);
  }
  printf("Solution x: ");
  fprint_dvector(stdout, n, x);

  /* solution check */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, mat_ptr(a), n, vec_ptr(x), 1,
              -1, vec_ptr(b), 1);
  norm = cblas_dnrm2(n, vec_ptr(b), 1);
  printf("|| A x - b || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: solution check\n");
    exit(255);
  }

  free_dmatrix(a);
  free_dvector(b);
  free_dmatrix(lu);
  free_dvector(x);
  free_ivector(ipiv);
  return 0;
}
