/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_H
#define ROKKO_SCALAPACK_H

#include <rokko/config.h>
#include <rokko/mangling.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SCALAPACK_pdelget ROKKO_GLOBAL(pdelget,PDELGET)
void SCALAPACK_pdelget(const char* scope, const char* top, double* alpha,
                       const double* A, const int* ia, const int* ja, const int* descA);

#define SCALAPACK_pdelset ROKKO_GLOBAL(pdelset,PDELSET)
void SCALAPACK_pdelset(double* A, const int* ia, const int* ja, const int* descA,
                       const double* alpha);

#define SCALAPACK_pdlamch ROKKO_GLOBAL(pdlamch,PDLAMCH)
double SCALAPACK_pdlamch(const int* icnt, const char* cmch);

#define SCALAPACK_pdlaprnt ROKKO_GLOBAL(pdlaprnt,PDLAPRNT)
void SCALAPACK_pdlaprnt(const int* m, const int* n, const double* A, const int* ia, const int* ja,
                        const int* descA, const int* irprnt, const int* icprnt, const char* cmatnm,
                        const int* nout, double* work);

#define SCALAPACK_pdsyev ROKKO_GLOBAL(pdsyev,PDSYEV)
void SCALAPACK_pdsyev(const char* jobz, const char* uplo, const int* n,
                      double* A, const int* ia, const int* ja, const int* descA,
                      double* w, double* Z, const int* iz, const int* jz, const int* descZ,
                      double* work, const int* lwork, int* info);

#define SCALAPACK_pdsyevd ROKKO_GLOBAL(pdsyevd,PDSYEVD)
void SCALAPACK_pdsyevd(const char* jobz, const char* uplo, const int* n,
                       double* A, const int* ia, const int* ja, const int* descA,
                       double* w,
                       double* Z, const int* iz, const int* jz, const int* descZ,
                       double* work, const int* lwork, int* iwork, const int* liwork, int* info);

#ifdef ROKKO_HAVE_PDSYEVR
#define SCALAPACK_pdsyevr ROKKO_GLOBAL(pdsyevr,PDSYEVR)
void SCALAPACK_pdsyevr(const char* jobz, const char* range, const char* uplo, const int* n,
                       double* A, const int* ia, const int* ja, const int* descA,
                       const double* vl, const double* vu, const int* il, const int* iu,
                       int* m, int* nz, double* w,
                       double* Z, const int* iz, const int* jz, const int* descZ,
                       double* work, const int* lwork, int* iwork, const int* liwork, int* info);
#endif

#define SCALAPACK_pdsyevx ROKKO_GLOBAL(pdsyevx,PDSYEVX)
void SCALAPACK_pdsyevx(const char* jobz, const char* range, const char* uplo, const int* n,
                       double* A, const int* iA, const int* jA, const int* descA,
                       const double* vl, const double* vu, const int* il, const int* iu,
                       const double* abstol, int* m, int* nZ, double* w, const double* orfac,
                       double* Z, const int* iZ, const int* jZ, const int* descZ,
                       double* work, const int* lwork, int* iwork, const int* liwork,
                       int* ifail, int* iclustr, double* gap, int* info);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_SCALAPACK_H
