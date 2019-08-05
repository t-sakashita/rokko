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
  
void cscalapack_pdelget(char scope, char top, double* alpha, const double* A, int ia, int ja,
                        const int* descA);

void cscalapack_pdelset(double* A, int ia, int ja, const int* descA, double alpha);

double cscalapack_pdlamch(int icnt, char cmch);

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

#ifdef ROKKO_HAVE_PDSYEVR
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

#ifdef ROKKO_HAVE_PDSYEVR
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

#ifdef __cplusplus
}
#endif

#endif // ROKKO_CSCALAPACK_H
