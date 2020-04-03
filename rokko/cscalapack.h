/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_CSCALAPACK_H
#define ROKKO_CSCALAPACK_H

#include <rokko/config.h>

#ifdef __cplusplus
extern "C" {
#endif

int cscalapack_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                        int ictxt, int lld);

int cscalapack_indxg2p(int indxglob, int nb, int iproc, int isrcproc, int nprocs);
  
int cscalapack_numroc(int n, int nb, int iproc, int isrcproc, int nprocs);

double cscalapack_pdelget(char scope, char top, const double* A, int ia, int ja, const int* descA);

void cscalapack_pdelset(double* A, int ia, int ja, const int* descA, double alpha);

double cscalapack_pdlamch(int icnt, char cmch);

double cscalapack_pdlange(char norm, int m, int n, const double* A, const int* descA);

double cscalapack_pdlange_work(char norm, int m, int n, const double* A, const int* descA,
                               double* work);

void cscalapack_pdlaprnt(int m, int n, const double* A, int ia, int ja, const int* descA,
                         int irprnt, int icprnt, const char* cmatnm, int nout, double* work);

int cscalapack_pdsyev_work(char jobz, char uplo, int n,
                           double* A, int ia, int ja, const int* descA,
                           double* w, double* Z, int iz, int jz, const int* descZ,
                           double* work, int lwork);

int cscalapack_pdsyevd_work(char jobz, char uplo, int n,
                            double* A, int ia, int ja, const int* descA,
                            double* w, double* Z, int iz, int jz, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork);

#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
int cscalapack_pdsyevr_work(char jobz, char range, char uplo, int n,
                            double* A, int ia, int ja, const int* descA,
                            double vl, double vu, int il, int iu,
                            int* m, int* nz, double* w,
                            double* Z, int iz, int jz, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork);
#endif

int cscalapack_pdsyevx_work(char jobz, char range, char uplo, int n,
                            double* A, int iA, int jA, const int* descA,
                            double vl, double vu, int il, int iu,
                            double abstol, int* m, int* nZ, double* w, double orfac,
                            double* Z, int iZ, int jZ, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork,
                            int* ifail, int* iclustr, double* gap);

int cscalapack_pdsyev(char jobz, char uplo, int n,
                      double* A, int ia, int ja, const int* descA,
                      double* w, double* Z, int iz, int jz, const int* descZ);

int cscalapack_pdsyevd(char jobz, char uplo, int n,
                       double* A, int ia, int ja, const int* descA,
                       double* w, double* Z, int iz, int jz, const int* descZ);

#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
int cscalapack_pdsyevr(char jobz, char range, char uplo, int n,
                       double* A, int ia, int ja, const int* descA,
                       double vl, double vu, int il, int iu,
                       int* m, int* nz, double* w,
                       double* Z, int iz, int jz, const int* descZ);
#endif

int cscalapack_pdsyevx(char jobz, char range, char uplo, int n,
                       double* A, int iA, int jA, const int* descA,
                       double vl, double vu, int il, int iu,
                       double abstol, int* m, int* nZ, double* w, double orfac,
                       double* Z, int iZ, int jZ, const int* descZ,
                       int* ifail, int* iclustr, double* gap);

float cscalapack_pselget(char scope, char top, const float* A, int ia, int ja, const int* descA);

void cscalapack_pselset(float* A, int ia, int ja, const int* descA, float alpha);

float cscalapack_pslamch(int icnt, char cmch);

float cscalapack_pslange(char norm, int m, int n, const float* A, const int* descA);

float cscalapack_pslange_work(char norm, int m, int n, const float* A, const int* descA,
                              float* work);

void cscalapack_pslaprnt(int m, int n, const float* A, int ia, int ja, const int* descA,
                         int irprnt, int icprnt, const char* cmatnm, int nout, float* work);

int cscalapack_pssyev_work(char jobz, char uplo, int n,
                           float* A, int ia, int ja, const int* descA,
                           float* w, float* Z, int iz, int jz, const int* descZ,
                           float* work, int lwork);

int cscalapack_pssyevd_work(char jobz, char uplo, int n,
                            float* A, int ia, int ja, const int* descA,
                            float* w, float* Z, int iz, int jz, const int* descZ,
                            float* work, int lwork, int* iwork, int liwork);

#ifdef ROKKO_HAVE_PSSYEVR
int cscalapack_pssyevr_work(char jobz, char range, char uplo, int n,
                            float* A, int ia, int ja, const int* descA,
                            float vl, float vu, int il, int iu,
                            int* m, int* nz, float* w,
                            float* Z, int iz, int jz, const int* descZ,
                            float* work, int lwork, int* iwork, int liwork);
#endif

int cscalapack_pssyevx_work(char jobz, char range, char uplo, int n,
                            float* A, int iA, int jA, const int* descA,
                            float vl, float vu, int il, int iu,
                            float abstol, int* m, int* nZ, float* w, float orfac,
                            float* Z, int iZ, int jZ, const int* descZ,
                            float* work, int lwork, int* iwork, int liwork,
                            int* ifail, int* iclustr, float* gap);

int cscalapack_pssyev(char jobz, char uplo, int n,
                      float* A, int ia, int ja, const int* descA,
                      float* w, float* Z, int iz, int jz, const int* descZ);

int cscalapack_pssyevd(char jobz, char uplo, int n,
                       float* A, int ia, int ja, const int* descA,
                       float* w, float* Z, int iz, int jz, const int* descZ);

#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
int cscalapack_pssyevr(char jobz, char range, char uplo, int n,
                       float* A, int ia, int ja, const int* descA,
                       float vl, float vu, int il, int iu,
                       int* m, int* nz, float* w,
                       float* Z, int iz, int jz, const int* descZ);
#endif

int cscalapack_pssyevx(char jobz, char range, char uplo, int n,
                       float* A, int iA, int jA, const int* descA,
                       float vl, float vu, int il, int iu,
                       float abstol, int* m, int* nZ, float* w, float orfac,
                       float* Z, int iZ, int jZ, const int* descZ,
                       int* ifail, int* iclustr, float* gap);
  
#ifdef __cplusplus
}
#endif

#endif // ROKKO_CSCALAPACK_H
