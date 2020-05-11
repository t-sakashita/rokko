/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_CSCALAPACK_H
#define ROKKO_CSCALAPACK_H

#include <rokko/config.h>
#include <lapacke.h>

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

#define CSCALAPACK_PSYEV_DECL(NAMES, NAMEL, TYPE) \
int cscalapack_ ## NAMES ## _work (char jobz, char uplo, int n, \
                                   TYPE* A, int ia, int ja, const int* descA, \
                                   TYPE* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                   TYPE* work, int lwork); \
int cscalapack_ ## NAMES (char jobz, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA, \
                          TYPE* w, TYPE* Z, int iz, int jz, const int* descZ);

#ifdef __cplusplus
extern "C" {
#endif

CSCALAPACK_PSYEV_DECL(pssyev, PSSYEV, float);
CSCALAPACK_PSYEV_DECL(pdsyev, PDSYEV, double);
CSCALAPACK_PSYEV_DECL(pcheev, PCHEEV, lapack_complex_float);
CSCALAPACK_PSYEV_DECL(pzheev, PZHEEV, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CSCALAPACK_PSYEV_DECL


#define CSCALAPACK_PSYEVD_DECL(NAMES, NAMEL, TYPE) \
int cscalapack_ ## NAMES ## _work (char jobz, char uplo, int n, \
                                   TYPE* A, int ia, int ja, const int* descA, \
                                   TYPE* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                   TYPE* work, int lwork, int* iwork, int liwork); \
int cscalapack_ ## NAMES (char jobz, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA, \
                          TYPE* w, TYPE* Z, int iz, int jz, const int* descZ);

#ifdef __cplusplus
extern "C" {
#endif

CSCALAPACK_PSYEVD_DECL(pssyevd, PSSYEVD, float);
CSCALAPACK_PSYEVD_DECL(pdsyevd, PDSYEVD, double);
CSCALAPACK_PSYEVD_DECL(pcheevd, PCHEEVD, lapack_complex_float);
CSCALAPACK_PSYEVD_DECL(pzheevd, PZHEEVD, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CSCALAPACK_PSYEVD_DECL


#define CSCALAPACK_PSYEVX_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_ ## NAMES ## _work(char jobz, char range, char uplo, int n, \
                                  TYPE* A, int ia, int ja, const int* descA, \
                                  TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                                  TYPE_REAL abstol, int* m, int* nZ, TYPE_REAL* w, TYPE_REAL orfac, \
                                  TYPE* Z, int iz, int jz, const int* descZ, \
                                  TYPE* work, int lwork, int* iwork, int liwork, \
                                  int* ifail, int* iclustr, TYPE_REAL* gap); \
int cscalapack_ ## NAMES (char jobz, char range, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA, \
                          TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                          TYPE_REAL abstol, int* m, int* nZ, TYPE_REAL* w, TYPE_REAL orfac, \
                          TYPE* Z, int iz, int jz, const int* descZ, \
                          int* ifail, int* iclustr, TYPE_REAL* gap);

#ifdef __cplusplus
extern "C" {
#endif

CSCALAPACK_PSYEVX_DECL(pssyevx, PSSYEVX, float, float);
CSCALAPACK_PSYEVX_DECL(pdsyevx, PDSYEVX, double, double);
CSCALAPACK_PSYEVX_DECL(pcheevx, PCHEEVX, lapack_complex_float, float);
CSCALAPACK_PSYEVX_DECL(pzheevx, PZHEEVX, lapack_complex_double, double);

#ifdef __cplusplus
}
#endif

#undef CSCALAPACK_PSYEVX_DECL


#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
#define CSCALAPACK_PSYEVR_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_## NAMES ##_work(char jobz, char range, char uplo, int n, \
                                TYPE* A, int ia, int ja, const int* descA, \
                                TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                                int* m, int* nz, \
                                TYPE_REAL* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                TYPE* work, int lwork, int* iwork, int liwork); \
int cscalapack_## NAMES (char jobz, char range, char uplo, int n, \
                         TYPE* A, int ia, int ja, const int* descA, \
                         TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                         int* m, int* nz, \
                         TYPE_REAL* w, TYPE* Z, int iz, int jz, const int* descZ);

#ifdef __cplusplus
extern "C" {
#endif

CSCALAPACK_PSYEVR_DECL(pssyevr, PSSYEVR, float, float);
CSCALAPACK_PSYEVR_DECL(pdsyevr, PDSYEVR, double, double);
CSCALAPACK_PSYEVR_DECL(pcheevr, PCHEEVR, lapack_complex_float, float);
CSCALAPACK_PSYEVR_DECL(pzheevr, PZHEEVR, lapack_complex_double, double);

#ifdef __cplusplus
}
#endif

#undef CSCALAPACK_PSYEVR_DECL
#endif

float cscalapack_pselget(char scope, char top, const float* A, int ia, int ja, const int* descA);

void cscalapack_pselset(float* A, int ia, int ja, const int* descA, float alpha);

float cscalapack_pslamch(int icnt, char cmch);

float cscalapack_pslange(char norm, int m, int n, const float* A, const int* descA);

float cscalapack_pslange_work(char norm, int m, int n, const float* A, const int* descA,
                              float* work);

void cscalapack_pslaprnt(int m, int n, const float* A, int ia, int ja, const int* descA,
                         int irprnt, int icprnt, const char* cmatnm, int nout, float* work);

int cscalapack_pdstebz(int ictxt, char range, char order, int n,
                       double vl, double vu, int il, int iu,
                       double abstol, const double* d, const double* e, int* m, int* nsplit,
                       double* w, int* iblock, int* isplit);

int cscalapack_pdstebz_work(int ictxt, char range, char order, int n,
                            double vl, double vu, int il, int iu,
                            double abstol, const double* d, const double* e, int* m, int* nsplit,
                            double* w, int* iblock, int* isplit,
                            double* work, int lwork, int* iwork, int liwork);

int cscalapack_pdstein(int n, const double* d, const double* e, int m,
                       double* w, const int* iblock, const int* isplit, double orfac,
                       double* Z, const int* iZ, const int* jZ, const int* descZ,
                       int* ifail, int* iclustr, double* gap);

int cscalapack_pdstein_work(int n, const double* d, const double* e, int m,
                            double* w, const int* iblock, const int* isplit, double orfac,
                            double* Z, const int* iZ, const int* jZ, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork,
                            int* ifail, int* iclustr, double* gap);

int cscalapack_pdstedc(char compz, int n, double* d, double* e,
                       double* Q, int iq, int jq, const int* descQ);

int cscalapack_pdstedc_work(char compz, int n, double* d, double* e,
                            double* Q, int iq, int jq, const int* descQ,
                            double* work, int lwork, int* iwork, int liwork);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_CSCALAPACK_H
