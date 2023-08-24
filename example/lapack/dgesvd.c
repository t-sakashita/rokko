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
#include <rokko/cmatrix.h>

int imin(int x, int y) { return (x < y) ? x : y; }
int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int m = 3;
  int n = 5;
  int r, lwork, i, j, info;
  double norm;
  double *s, *work;
  double **a, **u, **vt, **t;

  if (argc > 2) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
  }
  r = imin(m, n);

  /* generate matrix */
  a = alloc_dmatrix(m, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      mat_elem(a, i, j) = n - imax(i, j);
    }
  }
  printf("Matrix A: ");
  fprint_dmatrix(stdout, m, n, a);

  /* singular value decomposition */
  u = alloc_dmatrix(m, r);
  vt = alloc_dmatrix(r, n);
  s = alloc_dvector(r);
  t = alloc_dmatrix(m, n);
  cblas_dcopy(m * n, mat_ptr(a), 1, mat_ptr(t), 1);
  lwork = imax(3 * r + imax(m, n), 5 * r);
  work = alloc_dvector(lwork);
  info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', m, n, mat_ptr(t), m,
                             vec_ptr(s), mat_ptr(u), m, mat_ptr(vt), r,
                             vec_ptr(work), lwork);
  if (info != 0) {
    fprintf(stderr, "Error: dgesvd fails\n");
    exit(255);
  }
  free_dmatrix(t);
  printf("Matrix U: ");
  fprint_dmatrix(stdout, m, r, u);
  printf("Singular values: ");
  fprint_dvector(stdout, r, s);
  printf("Matrix Vt: ");
  fprint_dmatrix(stdout, r, n, vt);

  /* orthogonality check */
  t = alloc_dmatrix(r, r);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, r, r, m, 1, mat_ptr(u), m,
              mat_ptr(u), m, 0, mat_ptr(t), r);
  for (i = 0; i < r; ++i) mat_elem(t, i, i) -= 1;
  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', r, r, mat_ptr(t), r);
  printf("|| U^t U - I || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, r, r, n, 1, mat_ptr(vt), r,
              mat_ptr(vt), r, 0, mat_ptr(t), r);
  for (i = 0; i < r; ++i) mat_elem(t, i, i) -= 1;
  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', r, r, mat_ptr(t), r);
  printf("|| V^t V - I || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: orthogonality check\n");
    exit(255);
  }
  free_dmatrix(t);

  /* solution check */
  t = alloc_dmatrix(m, r);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, r, n, 1, mat_ptr(a), m,
              mat_ptr(vt), r, 0, mat_ptr(t), m);
  free_dmatrix(vt);
  vt = alloc_dmatrix(r, r);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, r, r, m, 1, mat_ptr(u), m,
              mat_ptr(t), m, 0, mat_ptr(vt), r);
  for (i = 0; i < r; ++i) mat_elem(vt, i, i) -= s[i];
  norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', r, r, mat_ptr(vt), r);
  printf("|| U^t A V - diag(S) || = %e\n", norm);
  if (norm > 1e-10) {
    fprintf(stderr, "Error: solution check\n");
    exit(255);
  }
  free_dmatrix(t);

  free_dmatrix(a);
  free_dmatrix(u);
  free_dmatrix(vt);
  free_dvector(s);
  return 0;
}
