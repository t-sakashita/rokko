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
  double **a, **ainv, **t;
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

  /* invert matrix */
  ainv = alloc_dmatrix(n, n);
  cblas_dcopy(n * n, MAT_PTR(a), 1, MAT_PTR(ainv), 1);
  ipiv = alloc_ivector(n);
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, MAT_PTR(ainv), n, VEC_PTR(ipiv));
  if (info != 0) {
    fprintf(stderr, "Error: dgetrf fails\n");
    exit(255);
  }
  info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, MAT_PTR(ainv), n, VEC_PTR(ipiv));
  if (info != 0) {
    fprintf(stderr, "Error: dgetri fails\n");
    exit(255);
  }
  printf("Matrix A^{-1}: ");
  fprint_dmatrix(stdout, n, n, ainv);

  /* solution check */
  t = alloc_dmatrix(n, n);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, MAT_PTR(ainv), n,
              MAT_PTR(a), n, 0, MAT_PTR(t), n);
  for (i = 0; i < n; ++i) MAT_ELEM(t, i, i) -= 1;
  norm2 = 0;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      norm2 = MAT_ELEM(t, i, j) * MAT_ELEM(t, i, j);
    }
  }
  printf("|| A^{-1} A - I ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: solution check\n");
    exit(255);
  }

  free_dmatrix(a);
  free_dmatrix(ainv);
  free_ivector(ipiv);
  free_dmatrix(t);
  return 0;
}
