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
  double *w;
  double complex **a, **u, **t;
  double complex alpha, beta;

  if (argc > 1) n = atoi(argv[1]);
    
  /* generate matrix */
  a = alloc_zmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < j; ++i) {
      mat_elem(a, i, j) = conj(mat_elem(a, j, i));
    }
    for (i = j; i < n; ++i) {
      mat_elem(a, i, j) = (n - imax(i, j)) + _Complex_I * (i - j);
    }
  }
  printf("Matrix A: ");
  fprint_zmatrix(stdout, n, n, a);

  /* diagonalization */
  u = alloc_zmatrix(n, n);
  cblas_zcopy(n * n, mat_ptr(a), 1, mat_ptr(u), 1);
  w = alloc_dvector(n);
  info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', n, mat_ptr(u), n, vec_ptr(w));
  if (info != 0) {
    fprintf(stderr, "Error: zheev fails with error code %d\n", info);
    exit(255);
  }
  printf("Eigenvalues: ");
  fprint_dvector(stdout, n, w);
  printf("Eigenvectors: ");
  fprint_zmatrix(stdout, n, n, u);

  /* orthogonality check */
  t = alloc_zmatrix(n, n);
  alpha = 1;
  beta = 0;
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, mat_ptr(u), n,
              mat_ptr(u), n, &beta, mat_ptr(t), n);
  for (i = 0; i < n; ++i) mat_elem(t, i, i) -= 1;
  norm = LAPACKE_zlange(LAPACK_COL_MAJOR, 'F', n, n, mat_ptr(t), n);
  printf("|| U^t U - I || = %e\n", norm);
  if (norm> 1e-10) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }

  /* eigenvalue check */
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, mat_ptr(a), n,
              mat_ptr(u), n, &beta, mat_ptr(t), n);
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, mat_ptr(u), n,
              mat_ptr(t), n, &beta, mat_ptr(a), n);
  for (i = 0; i < n; ++i) mat_elem(a, i, i) -= w[i];
  norm = LAPACKE_zlange(LAPACK_COL_MAJOR, 'F', n, n, mat_ptr(a), n);
  printf("|| U^t A U - diag(w) || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: eigenvalue check\n");
    exit(255);
  }

  free_zmatrix(a);
  free_zmatrix(u);
  free_dvector(w);
  free_zmatrix(t);
  return 0;
}
