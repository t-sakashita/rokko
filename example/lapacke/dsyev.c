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
  double *w;
  double **a, **u, **t;

  if (argc > 1) n = atoi(argv[1]);
    
  /* generate matrix */
  a = alloc_dmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      MAT_ELEM(a, i, j) = n - imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix(stdout, n, n, a);

  /* diagonalization */
  u = alloc_dmatrix(n, n);
  cblas_dcopy(n * n, MAT_PTR(a), 1, MAT_PTR(u), 1);
  w = alloc_dvector(n);
  info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, MAT_PTR(u), n, VEC_PTR(w));
  printf("Eigenvalues: ");
  fprint_dvector(stdout, n, w);
  printf("Eigenvectors: ");
  fprint_dmatrix(stdout, n, n, u);

  /* orthogonality check */
  t = alloc_dmatrix(n, n);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, n, 1, MAT_PTR(u), n,
              MAT_PTR(u), n, 0, MAT_PTR(t), n);
  for (i = 0; i < n; ++i) MAT_ELEM(t, i, i) -= 1;
  norm2 = 0;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      norm2 = MAT_ELEM(t, i, j) * MAT_ELEM(t, i, j);
    }
  }
  printf("|| U^t U - I ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }

  /* eigenvalue check */
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, MAT_PTR(a), n,
              MAT_PTR(u), n, 0, MAT_PTR(t), n);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, n, 1, MAT_PTR(u), n,
              MAT_PTR(t), n, 0, MAT_PTR(a), n);
  for (i = 0; i < n; ++i) MAT_ELEM(a, i, i) -= w[i];
  norm2 = 0;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      norm2 = MAT_ELEM(a, i, j) * MAT_ELEM(a, i, j);
    }
  }
  printf("|| U^t A U - diag(w) ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: eigenvalue check\n");
    exit(255);
  }

  free_dmatrix(a);
  free_dmatrix(u);
  free_dvector(w);
  free_dmatrix(t);
  return 0;
}
