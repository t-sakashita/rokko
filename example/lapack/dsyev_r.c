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

/* dsyev using row-major matrix */

#include <cblas.h>
#include <lapacke.h>
#include <cmatrix.h>

int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int n = 5;
  int i, j, info;
  double norm;
  double *w;
  double **a, **u, **t;

  if (argc > 1) n = atoi(argv[1]);

  /* generate matrix */
  a = alloc_dmatrix_r(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      a[i][j] = n - imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix_r(stdout, n, n, a);

  /* diagonalization */
  u = alloc_dmatrix_r(n, n);
  cblas_dcopy(n * n, mat_ptr(a), 1, mat_ptr(u), 1);
  w = alloc_dvector(n);
  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, mat_ptr(u), n, vec_ptr(w));
  if (info != 0) {
    fprintf(stderr, "Error: dsyev fails with error code %d\n", info);
    exit(255);
  }
  printf("Eigenvalues: ");
  fprint_dvector(stdout, n, w);
  printf("Eigenvectors: ");
  fprint_dmatrix_r(stdout, n, n, u);

  /* orthogonality check */
  t = alloc_dmatrix_r(n, n);
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1, mat_ptr(u), n,
              mat_ptr(u), n, 0, mat_ptr(t), n);
  for (i = 0; i < n; ++i) t[i][i] -= 1;
  norm = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'F', n, n, mat_ptr(t), n);
  printf("|| U^t U - I || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }

  /* eigenvalue check */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, mat_ptr(a), n,
              mat_ptr(u), n, 0, mat_ptr(t), n);
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1, mat_ptr(u), n,
              mat_ptr(t), n, 0, mat_ptr(a), n);
  for (i = 0; i < n; ++i) a[i][i] -= w[i];
  norm = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'F', n, n, mat_ptr(a), n);
  printf("|| U^t A U - diag(w) || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: eigenvalue check\n");
    exit(255);
  }

  free_dmatrix(a);
  free_dmatrix(u);
  free_dvector(w);
  free_dmatrix(t);
  return 0;
}
