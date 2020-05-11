/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_SCALAPACK_INTERFACE_H
#define ROKKO_SCALAPACK_SCALAPACK_INTERFACE_H

#include <rokko/config.h>
#include <rokko/mangling.h>
#include <lapacke.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SCALAPACK_descinit ROKKO_GLOBAL(descinit,DESCINIT)
void SCALAPACK_descinit(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                        const int* irsrc, const int* icsrc, const int* ictxt, const int* lld,
                        int* info);
  
#define SCALAPACK_indxg2p ROKKO_GLOBAL(indxg2p,INDXG2P)
int SCALAPACK_indxg2p(int* indxglob, int* nb, int* iproc, int* isrcproc, int* nprocs);

#define SCALAPACK_numroc ROKKO_GLOBAL(numroc,NUMROC)
int SCALAPACK_numroc(int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);

#define SCALAPACK_pdelget ROKKO_GLOBAL(pdelget,PDELGET)
void SCALAPACK_pdelget(const char* scope, const char* top, double* alpha,
                       const double* A, const int* ia, const int* ja, const int* descA);

#define SCALAPACK_pdelset ROKKO_GLOBAL(pdelset,PDELSET)
void SCALAPACK_pdelset(double* A, const int* ia, const int* ja, const int* descA,
                       const double* alpha);

#define SCALAPACK_pdlamch ROKKO_GLOBAL(pdlamch,PDLAMCH)
double SCALAPACK_pdlamch(const int* icnt, const char* cmch);

#define SCALAPACK_pdlange ROKKO_GLOBAL(pdlange,PDLANGE)
double SCALAPACK_pdlange(const char* norm, int* m, int* n, const double* A,
                         const int* ia, const int* ja, const int* descA, double* work);

#define SCALAPACK_pdlaprnt ROKKO_GLOBAL(pdlaprnt,PDLAPRNT)
void SCALAPACK_pdlaprnt(const int* m, const int* n, const double* A, const int* ia, const int* ja,
                        const int* descA, const int* irprnt, const int* icprnt, const char* cmatnm,
                        const int* nout, double* work);

#define SCALAPACK_PSYEV_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* uplo, const int* n, \
                                 TYPE* A, const int* ia, const int* ja, const int* descA, \
                                 TYPE* w, TYPE* Z, const int* iz, const int* jz, const int* descZ, \
                                 TYPE* work, const int* lwork, int* info);

#ifdef __cplusplus
extern "C" {
#endif

SCALAPACK_PSYEV_DECL(pssyev, PSSYEV, float);
SCALAPACK_PSYEV_DECL(pdsyev, PDSYEV, double);
SCALAPACK_PSYEV_DECL(pcheev, PCHEEV, lapack_complex_float);
SCALAPACK_PSYEV_DECL(pzheev, PZHEEV, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef SCALAPACK_PSYEV_DECL


#define SCALAPACK_PSYEVD_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* uplo, const int* n, \
  TYPE* A, const int* ia, const int* ja, const int* descA, \
  TYPE* w, \
  TYPE* Z, const int* iz, const int* jz, const int* descZ, \
  TYPE* work, const int* lwork, int* iwork, const int* liwork, int* info);

#ifdef __cplusplus
extern "C" {
#endif

SCALAPACK_PSYEVD_DECL(pssyevd, PSSYEVD, float);
SCALAPACK_PSYEVD_DECL(pdsyevd, PDSYEVD, double);
SCALAPACK_PSYEVD_DECL(pcheevd, PCHEEVD, lapack_complex_float);
SCALAPACK_PSYEVD_DECL(pzheevd, PZHEEVD, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef SCALAPACK_PSYEVD_DECL


#define SCALAPACK_PSYEVX_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* range, const char* uplo, const int* n, \
  TYPE* A, const int* iA, const int* jA, const int* descA, \
  const TYPE_REAL* vl, const TYPE_REAL* vu, const int* il, const int* iu, \
  const TYPE_REAL* abstol, int* m, int* nZ, TYPE_REAL* w, const TYPE_REAL* orfac, \
  TYPE* Z, const int* iZ, const int* jZ, const int* descZ, \
  TYPE* work, const int* lwork, int* iwork, const int* liwork, \
  int* ifail, int* iclustr, TYPE_REAL* gap, int* info);

#ifdef __cplusplus
extern "C" {
#endif

SCALAPACK_PSYEVX_DECL(pssyevx, PSSYEVX, float, float);
SCALAPACK_PSYEVX_DECL(pdsyevx, PDSYEVX, double, double);
SCALAPACK_PSYEVX_DECL(pcheevx, PCHEEVX, lapack_complex_float, float);
SCALAPACK_PSYEVX_DECL(pzheevx, PZHEEVX, lapack_complex_double, double);

#ifdef __cplusplus
}
#endif

#undef SCALAPACK_PSYEVX_DECL


#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
#define SCALAPACK_pdsyevr ROKKO_GLOBAL(pdsyevr,PDSYEVR)
void SCALAPACK_pdsyevr(const char* jobz, const char* range, const char* uplo, const int* n,
                       double* A, const int* ia, const int* ja, const int* descA,
                       const double* vl, const double* vu, const int* il, const int* iu,
                       int* m, int* nz, double* w,
                       double* Z, const int* iz, const int* jz, const int* descZ,
                       double* work, const int* lwork, int* iwork, const int* liwork, int* info);
#endif

#define SCALAPACK_pselget ROKKO_GLOBAL(pselget,PSELGET)
void SCALAPACK_pselget(const char* scope, const char* top, float* alpha,
                       const float* A, const int* ia, const int* ja, const int* descA);

#define SCALAPACK_pselset ROKKO_GLOBAL(pselset,PSELSET)
void SCALAPACK_pselset(float* A, const int* ia, const int* ja, const int* descA,
                       const float* alpha);

#define SCALAPACK_pslamch ROKKO_GLOBAL(pslamch,PSLAMCH)
float SCALAPACK_pslamch(const int* icnt, const char* cmch);

#define SCALAPACK_pslange ROKKO_GLOBAL(pslange,PSLANGE)
float SCALAPACK_pslange(const char* norm, int* m, int* n, const float* A,
                        const int* ia, const int* ja, const int* descA, float* work);

#define SCALAPACK_pslaprnt ROKKO_GLOBAL(pslaprnt,PSLAPRNT)
void SCALAPACK_pslaprnt(const int* m, const int* n, const float* A, const int* ia, const int* ja,
                        const int* descA, const int* irprnt, const int* icprnt, const char* cmatnm,
                        const int* nout, float* work);

#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
#define SCALAPACK_pssyevr ROKKO_GLOBAL(pssyevr,PSSYEVR)
void SCALAPACK_pssyevr(const char* jobz, const char* range, const char* uplo, const int* n,
                       float* A, const int* ia, const int* ja, const int* descA,
                       const float* vl, const float* vu, const int* il, const int* iu,
                       int* m, int* nz, float* w,
                       float* Z, const int* iz, const int* jz, const int* descZ,
                       float* work, const int* lwork, int* iwork, const int* liwork, int* info);
#endif

#define SCALAPACK_pdstebz ROKKO_GLOBAL(pdstebz,PDSTEBZ)
void SCALAPACK_pdstebz(const int* ictxt, const char* range, const char* order, const int* n,
                       const double* vl, const double* vu, const int* il, const int* iu,
                       const double* abstol, const double* d, const double* e, int* m, int* nsplit,
                       double* w, int* iblock, int* isplit,
                       double* work, const int* lwork, int* iwork, const int* liwork, int* info);

#define SCALAPACK_pdstein ROKKO_GLOBAL(pdstein,PDSTEIN)
void SCALAPACK_pdstein(const int* n, const double* d, const double* e, const int* m,
                       double* w, const int* iblock, const int* isplit, const double* orfac,
                       double* Z, const int* iZ, const int* jZ, const int* descZ,
                       double* work, const int* lwork, int* iwork, const int* liwork,
                       int* ifail, int* iclustr, double* gap, int* info);

#define SCALAPACK_pdstedc ROKKO_GLOBAL(pdstedc,PDSTEDC)
void SCALAPACK_pdstedc(const char* compz, const int* n, double* d, double* e,
                       double* Q, const int* iq, const int* jq, const int* descQ,
                       double* work, const int* lwork, int* iwork, const int* liwork, int* info);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_SCALAPACK_SCALAPACK_INTERFACE_H
