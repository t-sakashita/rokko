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

int imin(int x, int y) { return (x < y) ? x : y; }
int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int m = 3;
  int n = 5;
  int r, lwork, i, j, info;
  double norm2;
  double *s, *work;
  double **a, **u, **vt, **t;

  if (argc > 2) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
  }
  r = imin(m, n);
  
  a = alloc_dmatrix(m, n);
  u = alloc_dmatrix(m, r);
  vt = alloc_dmatrix(r, n);
  s = alloc_dvector(r);

  /* generate matrix */
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      MAT_ELEM(a, i, j) = n - imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix(stdout, m, n, a);

  /* singular value decomposition */
  lwork = imax(3 * r + imax(m, n), 5 * r);
  work = alloc_dvector(lwork);
  t = alloc_dmatrix(m, n);
  cblas_dcopy(m * n, MAT_PTR(a), 1, MAT_PTR(t), 1);
  info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S',
                             m, n, MAT_PTR(t), m,
                             VEC_PTR(s), MAT_PTR(u), m, MAT_PTR(vt), r,
                             VEC_PTR(work), lwork);
  free_dmatrix(t);
  printf("Matrix U: ");
  fprint_dmatrix(stdout, m, r, u);
  printf("Singular values: ");
  fprint_dvector(stdout, r, s);
  printf("Matrix Vt: ");
  fprint_dmatrix(stdout, r, n, vt);

  /* orthogonality check */
  t = alloc_dmatrix(r, r);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, r, r, m,
              1, MAT_PTR(u), m, MAT_PTR(u), m, 0, MAT_PTR(t), r);
  for (i = 0; i < r; ++i) MAT_ELEM(t, i, i) -= 1;
  norm2 = 0;
  for (j = 0; j < r; ++j) {
    for (i = 0; i < r; ++i) {
      norm2 = MAT_ELEM(t, i, j) * MAT_ELEM(t, i, j);
    }
  }
  printf("|| U^t U - I ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, r, r, n,
              1, MAT_PTR(vt), r, MAT_PTR(vt), r, 0, MAT_PTR(t), r);
  for (i = 0; i < r; ++i) MAT_ELEM(t, i, i) -= 1;
  norm2 = 0;
  for (j = 0; j < r; ++j) {
    for (i = 0; i < r; ++i) {
      norm2 = MAT_ELEM(t, i, j) * MAT_ELEM(t, i, j);
    }
  }
  printf("|| V^t V - I ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }
  free_dmatrix(t);

  /* solution check */
  t = alloc_dmatrix(m, r);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, r, n,
              1, MAT_PTR(a), m, MAT_PTR(vt), r, 0, MAT_PTR(t), m);
  free_dmatrix(vt);
  vt = alloc_dmatrix(r, r);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, r, r, m,
              1, MAT_PTR(u), m, MAT_PTR(t), m, 0, MAT_PTR(vt), r);
  for (i = 0; i < r; ++i) MAT_ELEM(vt, i, i) -= s[i];
  norm2 = 0;
  for (j = 0; j < r; ++j) {
    for (i = 0; i < r; ++i) {
      norm2 = MAT_ELEM(vt, i, j) * MAT_ELEM(vt, i, j);
    }
  }
  printf("|| U^t A V - diag(S) ||^2 = %e\n", norm2);
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: solution check\n");
    exit(255);
  }

  free_dvector(s);
  free_dvector(work);
  free_dmatrix(a);
  free_dmatrix(u);
  free_dmatrix(vt);
  free_dmatrix(t);
  return 0;
}
