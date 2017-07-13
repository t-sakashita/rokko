/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef CMATRIX_H
#define CMATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#undef I

#define vec_ptr(vec) &(vec)[0]
#define mat_ptr(mat) &(mat)[0][0]
#define mat_elem(mat, i, j) (mat)[j][i]

#ifdef __cplusplus
extern "C" {
#endif

/*********
  vectors
**********/

/* allocate vector of int */
static inline int *alloc_ivector(int n) {
  int *vec;
  vec = (int*)calloc(n, sizeof(int));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_ivector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of float */
static inline float *alloc_svector(int n) {
  float *vec;
  vec = (float*)calloc(n, sizeof(float));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_svector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of double */
static inline double *alloc_dvector(int n) {
  double *vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of float complex */
static inline float complex *alloc_cvector(int n) {
  float complex *vec;
  vec = (float complex*)calloc(n, sizeof(float complex));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_cvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of double complex */
static inline double complex *alloc_zvector(int n) {
  double complex *vec;
  vec = (double complex*)calloc(n, sizeof(double complex));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_zvector\n");
    exit(1);
  }
  return vec;
}

/* deallocate vector of int */
static inline void free_ivector(int *vec) {
  free(vec);
  vec = NULL;
}

/* deallocate vector of float */
static inline void free_svector(float *vec) {
  free(vec);
  vec = NULL;
}

/* deallocate vector of double */
static inline void free_dvector(double *vec) {
  free(vec);
  vec = NULL;
}

/* deallocate vector of float complex */
static inline void free_cvector(float complex *vec) {
  free(vec);
  vec = NULL;
}

/* deallocate vector of double complex */
static inline void free_zvector(double complex *vec) {
  free(vec);
  vec = NULL;
}

/* print out vector of int */
static inline void fprintf_ivector(FILE *fp, const char* fmt, int n, int *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, fmt, vec[i]);
  fprintf(fp, "\n");
}
static inline void fprint_ivector(FILE *fp, int n, int *vec) {
  fprintf_ivector(fp, "%d ", n, vec);
}

/* print out vector of float */
static inline void fprintf_svector(FILE *fp, const char* fmt, int n, float *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, fmt, vec[i]);
  fprintf(fp, "\n");
}
static inline void fprint_svector(FILE *fp, int n, float *vec) {
  fpintf_svector(fp, "%e ", n, vec);
}

/* print out vector of double */
static inline void fprintf_dvector(FILE *fp, const char* fmt, int n, double *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, fmt, vec[i]);
  fprintf(fp, "\n");
}
static inline void fprint_dvector(FILE *fp, int n, double *vec) {
  fprintf_dvector(fp, "%20.14e ", n, vec);
}

/* print out vector of float complex */
static inline void fprintf_cvector(FILE *fp, const char* fmt, int n, float complex *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, fmt, crealf(vec[i]), cimagf(vec[i]));
  fprintf(fp, "\n");
}
static inline void fprint_cvector(FILE *fp, int n, float complex *vec) {
  fpintf_cvector(fp, "(%e,%e) ", n, vec);
}

/* print out vector of double complex */
static inline void fprintf_zvector(FILE *fp, const char* fmt, int n, double complex *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, fmt, creal(vec[i]), cimag(vec[i]));
  fprintf(fp, "\n");
}
static inline void fprint_zvector(FILE *fp, int n, double complex *vec) {
  fprintf_zvector(fp, "(%20.14e,%20.14e) ", n, vec);
}

/* read vector of int from file */
static inline void read_ivector(FILE *fp, int *n, int **vec) {
  int i;
  fscanf(fp, "%d", n);
  *vec = alloc_ivector(*n);
  for (i = 0; i < *n; ++i) fscanf(fp, "%d", &(*vec)[i]);
}

/* read vector of float from file */
static inline void read_svector(FILE *fp, int *n, float **vec) {
  int i;
  fscanf(fp, "%d", n);
  *vec = alloc_svector(*n);
  for (i = 0; i < *n; ++i) fscanf(fp, "%e", &(*vec)[i]);
}

/* read vector of double from file */
static inline void read_dvector(FILE *fp, int *n, double **vec) {
  int i;
  fscanf(fp, "%d", n);
  *vec = alloc_dvector(*n);
  for (i = 0; i < *n; ++i) fscanf(fp, "%le", &(*vec)[i]);
}

/* read vector of float complex from file */
static inline void read_cvector(FILE *fp, int *n, float complex **vec) {
  int i;
  float v_re, v_im;
  fscanf(fp, "%d", n);
  *vec = alloc_cvector(*n);
  for (i = 0; i < *n; ++i) {
    fscanf(fp, "(%e,%e)%*c", &v_re, &v_im);
    (*vec)[i] = v_re + _Complex_I * v_im;
  }
}

/* read vector of double from file */
static inline void read_zvector(FILE *fp, int *n, double complex **vec) {
  int i;
  double v_re, v_im;
  fscanf(fp, "%d", n);
  *vec = alloc_zvector(*n);
  for (i = 0; i < *n; ++i) {
    fscanf(fp, "(%le,%le)%*c", &v_re, &v_im);
    (*vec)[i] = v_re + _Complex_I * v_im;
  }
}

/**********
  matrices
***********/

/* allocate m x n column-major matrix of float */
static inline float **alloc_smatrix(int m, int n) {
  int i;
  float **mat;
  mat = (float**)malloc((size_t)(n * sizeof(float*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_smatrix\n");
    exit(1);
  }
  mat[0] = (float*)malloc((size_t)(m * n * sizeof(float)));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_smatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n row-major matrix of float */
static inline float **alloc_smatrix_r(int m, int n) { return alloc_smatrix(n, m); }

/* allocate m x n column-major matrix of double */
static inline double **alloc_dmatrix(int m, int n) {
  int i;
  double **mat;
  mat = (double**)malloc((size_t)(n * sizeof(double*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  mat[0] = (double*)calloc(m * n, sizeof(double));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n row-major matrix of double */
static inline double **alloc_dmatrix_r(int m, int n) { return alloc_dmatrix(n, m); }

/* allocate m x n column-major matrix of float complex */
static inline float complex **alloc_cmatrix(int m, int n) {
  int i;
  float complex **mat;
  mat = (float complex**)malloc((size_t)(n * sizeof(float complex*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_cmatrix\n");
    exit(1);
  }
  mat[0] = (float complex*)malloc((size_t)(m * n * sizeof(float complex)));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_cmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n row-major matrix of float complex */
static inline float complex **alloc_cmatrix_r(int m, int n) { return alloc_cmatrix(n, m); }

/* allocate m x n column-major matrix of double complex */
static inline double complex **alloc_zmatrix(int m, int n) {
  int i;
  double complex **mat;
  mat = (double complex**)malloc((size_t)(n * sizeof(double complex*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_zmatrix\n");
    exit(1);
  }
  mat[0] = (double complex*)calloc(m * n, sizeof(double complex));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_zmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n row-major matrix of double complex */
static inline double complex **alloc_zmatrix_r(int m, int n) { return alloc_zmatrix(n, m); }

/* deallocate matrix of float */
static inline void free_smatrix(float **mat) {
  free(mat[0]);
  free(mat);
  mat = NULL;
}

/* deallocate matrix of double */
static inline void free_dmatrix(double **mat) {
  free(mat[0]);
  free(mat);
  mat = NULL;
}

/* deallocate matrix of float complex */
static inline void free_cmatrix(float complex **mat) {
  free(mat[0]);
  free(mat);
  mat = NULL;
}

/* deallocate matrix of double complex */
static inline void free_zmatrix(double complex **mat) {
  free(mat[0]);
  free(mat);
  mat = NULL;
}

/* print out column-major matrix of float */
static inline void fprint_smatrixf(FILE *fp, const char *fmt, int m, int n, float **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, fmt, mat_elem(mat, i, j));
    fprintf(fp, "\n");
  }
}
static inline void fprint_smatrix(FILE *fp, int m, int n, float **mat) {
  fprint_smatrixf(fp, "%e ", m, n, mat);
}

/* print out row-major matrix of float */
static inline void fprint_smatrixf_r(FILE *fp, const char *fmt, int m, int n, float **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, fmt, mat[i][j]);
    fprintf(fp, "\n");
  }
}
static inline void fprint_smatrix_r(FILE *fp, int m, int n, float **mat) {
  fprint_smatrixf_r(fp, "%e ", m, n, mat);
}

/* print out column-major matrix of double */
static inline void fprint_dmatrixf(FILE *fp, const char *fmt, int m, int n, double **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, fmt, mat_elem(mat, i, j));
    fprintf(fp, "\n");
  }
}
static inline void fprint_dmatrix(FILE *fp, int m, int n, double **mat) {
  fprint_dmatrixf(fp, "%20.14e ", m, n, mat);
}
/* print out row-major matrix of double */
static inline void fprint_dmatrixf_r(FILE *fp, const char *fmt, int m, int n, double **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, fmt, mat[i][j]);
    fprintf(fp, "\n");
  }
}
static inline void fprint_dmatrix_r(FILE *fp, int m, int n, double **mat) {
  fprint_dmatrixf_r(fp, "%20.14e ", m, n, mat);
}

/* print out column-major matrix of float complex */
static inline void fprint_cmatrixf(FILE *fp, const char *fmt, int m, int n, float complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j)
      fprintf(fp, fmt, crealf(mat_elem(mat, i, j)), cimagf(mat_elem(mat, i, j)));
    fprintf(fp, "\n");
  }
}
static inline void fprint_cmatrix(FILE *fp, int m, int n, float complex **mat) {
  fprint_cmatrixf(fp, "(%e,%e) ", m, n, mat);
}

/* print out row-major matrix of float complex */
static inline void fprint_cmatrixf_r(FILE *fp, const char *fmt, int m, int n, float complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j)
      fprintf(fp, fmt, crealf(mat[i][j]), cimagf(mat[i][j]));
    fprintf(fp, "\n");
  }
}
static inline void fprint_cmatrix_r(FILE *fp, int m, int n, float complex **mat) {
  fprint_cmatrixf_r(fp, "(%e,%e) ", m, n, mat);
}

/* print out column-major matrix of double complex */
static inline void fprint_zmatrixf(FILE *fp, const char *fmt, int m, int n, double complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j)
      fprintf(fp, fmt, creal(mat_elem(mat, i, j)), cimag(mat_elem(mat, i, j)));
    fprintf(fp, "\n");
  }
}
static inline void fprint_zmatrix(FILE *fp, int m, int n, double complex **mat) {
  fprint_zmatrixf(fp, "(%20.14e,%20.14e) ", m, n, mat);
}

/* print out row-major matrix of double complex */
static inline void fprint_zmatrixf_r(FILE *fp, const char *fmt, int m, int n, double complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j)
      fprintf(fp, fmt, creal(mat[i][j]), cimag(mat[i][j]));
    fprintf(fp, "\n");
  }
}
static inline void fprint_zmatrix_r(FILE *fp, int m, int n, double complex **mat) {
  fprint_zmatrixf_r(fp, "(%20.14e,%20.14e) ", m, n, mat);
}

/* read column-major matrix of float from file */
static inline void read_smatrix(FILE *fp, int *m, int *n, float ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_smatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%f", &mat_elem(*mat, i, j));
  }
}

/* read row-major matrix of float from file */
static inline void read_smatrix_r(FILE *fp, int *m, int *n, float ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_smatrix_r(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%f", &(*mat)[i][j]);
  }
}

/* read column-major matrix of double from file */
static inline void read_dmatrix(FILE *fp, int *m, int *n, double ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_dmatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%lf", &mat_elem(*mat, i, j));
  }
}

/* read row-major matrix of double from file */
static inline void read_dmatrix_r(FILE *fp, int *m, int *n, double ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_dmatrix_r(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%lf", &(*mat)[i][j]);
  }
}

/* read column-major matrix of float complex from file */
static inline void read_cmatrix(FILE *fp, int *m, int *n, float complex ***mat) {
  int i, j;
  float v_re, v_im;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_cmatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) {
      fscanf(fp, "(%f,%f)%*c", v_re, v_im);
      mat_elem(*mat, i, j) = v_re + _Complex_I * v_im;
    }
  }
}

/* read row-major matrix of float complex from file */
static inline void read_cmatrix_r(FILE *fp, int *m, int *n, float complex ***mat) {
  int i, j;
  float v_re, v_im;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_cmatrix_r(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) {
      fscanf(fp, "(%f,%f)%*c", v_re, v_im);
      (*mat)[i][j] = v_re + _Complex_I * v_im;
    }
  }
}

/* read column-major matrix of double complex from file */
static inline void read_zmatrix(FILE *fp, int *m, int *n, double complex ***mat) {
  int i, j;
  double v_re, v_im;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_zmatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) {
      fscanf(fp, "(%lf,%lf)%*c", v_re, v_im);
      mat_elem(*mat, i, j) = v_re + _Complex_I * v_im;
    }
  }
}

/* read row-major matrix of double complex from file */
static inline void read_zmatrix_r(FILE *fp, int *m, int *n, double complex ***mat) {
  int i, j;
  double v_re, v_im;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_zmatrix_r(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) {
      fscanf(fp, "(%lf,%lf)%*c", v_re, v_im);
      (*mat)[i][j] = v_re + _Complex_I * v_im;
    }
  }
}

#ifdef __cplusplus
}
#endif

#endif
