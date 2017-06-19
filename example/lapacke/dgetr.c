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
  double norm2;
  double *x, *b;
  double **a, **lu;
  int *ipiv;

  if (argc > 1) n = atoi(argv[1]);

  // generate matrix and rhs vector
  a = alloc_dmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      MAT_ELEM(a, i, j) = n - imax(i, j);
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
  cblas_dcopy(n * n, MAT_PTR(a), 1, MAT_PTR(lu), 1);
  x = alloc_dvector(n);
  cblas_dcopy(n, VEC_PTR(b), 1, VEC_PTR(x), 1);
  ipiv = alloc_ivector(n);
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, MAT_PTR(lu), n, VEC_PTR(ipiv));
  if (info != 0) {
    fprintf(stderr, "Error: dgetrf fails\n");
    exit(255);
  }
  info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'n', n, 1, MAT_PTR(lu), n, VEC_PTR(ipiv),
                        VEC_PTR(x), n);
  if (info != 0) {
    fprintf(stderr, "Error: dgetrs fails\n");
    exit(255);
  }
  printf("Solution x: ");
  fprint_dvector(stdout, n, x);

  /* solution check */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, MAT_PTR(a), n, VEC_PTR(x), 1,
              -1, VEC_PTR(b), 1);
  norm2 = 0;
  for (i = 0; i < n; ++i) norm2 = b[i] * b[i];
  printf("|| A x - b ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
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
