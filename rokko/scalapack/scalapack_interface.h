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

#define SCALAPACK_PSYEV_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* uplo, const int* n, \
                                 TYPE* A, const int* ia, const int* ja, const int* descA, \
                                 TYPE_REAL* w, TYPE* Z, const int* iz, const int* jz, const int* descZ, \
                                 TYPE* work, const int* lwork, int* info);

SCALAPACK_PSYEV_DECL(pssyev, PSSYEV, float, float);
SCALAPACK_PSYEV_DECL(pdsyev, PDSYEV, double, double);
SCALAPACK_PSYEV_DECL(pcheev, PCHEEV, lapack_complex_float, float);
SCALAPACK_PSYEV_DECL(pzheev, PZHEEV, lapack_complex_double, double);

#undef SCALAPACK_PSYEV_DECL


#define SCALAPACK_PSYEVD_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* uplo, const int* n, \
  TYPE* A, const int* ia, const int* ja, const int* descA, \
  TYPE_REAL* w, \
  TYPE* Z, const int* iz, const int* jz, const int* descZ, \
  TYPE* work, const int* lwork, int* iwork, const int* liwork, int* info);

SCALAPACK_PSYEVD_DECL(pssyevd, PSSYEVD, float, float);
SCALAPACK_PSYEVD_DECL(pdsyevd, PDSYEVD, double, double);
SCALAPACK_PSYEVD_DECL(pcheevd, PCHEEVD, lapack_complex_float, float);
SCALAPACK_PSYEVD_DECL(pzheevd, PZHEEVD, lapack_complex_double, double);

#undef SCALAPACK_PSYEVD_DECL


#define SCALAPACK_PSYEVX_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* range, const char* uplo, const int* n, \
  TYPE* A, const int* iA, const int* jA, const int* descA, \
  const TYPE_REAL* vl, const TYPE_REAL* vu, const int* il, const int* iu, \
  const TYPE_REAL* abstol, int* m, int* nZ, TYPE_REAL* w, const TYPE_REAL* orfac, \
  TYPE* Z, const int* iZ, const int* jZ, const int* descZ, \
  TYPE* work, const int* lwork, int* iwork, const int* liwork, \
  int* ifail, int* iclustr, TYPE_REAL* gap, int* info);

SCALAPACK_PSYEVX_DECL(pssyevx, PSSYEVX, float, float);
SCALAPACK_PSYEVX_DECL(pdsyevx, PDSYEVX, double, double);
SCALAPACK_PSYEVX_DECL(pcheevx, PCHEEVX, lapack_complex_float, float);
SCALAPACK_PSYEVX_DECL(pzheevx, PZHEEVX, lapack_complex_double, double);

#undef SCALAPACK_PSYEVX_DECL

#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
#define SCALAPACK_PSYEVR_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* jobz, const char* range, const char* uplo, const int* n, \
  TYPE* A, const int* iA, const int* jA, const int* descA, \
  const TYPE_REAL* vl, const TYPE_REAL* vu, const int* il, const int* iu, \
  int* m, int* nz, \
  TYPE_REAL* w, TYPE* Z, const int* iz, const int* jz, const int* descZ, \
  TYPE* work, const int* lwork, int* iwork, const int* liwork, int* info);

SCALAPACK_PSYEVR_DECL(pssyevr, PSSYEVR, float, float);
SCALAPACK_PSYEVR_DECL(pdsyevr, PDSYEVR, double, double);
SCALAPACK_PSYEVR_DECL(pcheevr, PCHEEVR, lapack_complex_float, float);
SCALAPACK_PSYEVR_DECL(pzheevr, PZHEEVR, lapack_complex_double, double);

#undef SCALAPACK_PSYEVR_DECL
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


#define SCALAPACK_PSTEBZ_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* ictxt, const char* range, const char* order, const int* n, \
                                 const TYPE* vl, const TYPE* vu, const int* il, const int* iu, \
                                 const TYPE* abstol, const TYPE* d, const TYPE* e, int* m, int* nsplit, \
                                 TYPE* w, int* iblock, int* isplit, \
                                 TYPE* work, const int* lwork, int* iwork, const int* liwork, int* info);

SCALAPACK_PSTEBZ_DECL(psstebz, PSSTEBZ, float);
SCALAPACK_PSTEBZ_DECL(pdstebz, PDSTEBZ, double);

#undef SCALAPACK_PSTEBZ_DECL


#define SCALAPACK_PSTEIN_DECL(NAMES, NAMEL, TYPE, TYPE_REAL) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* n, const TYPE_REAL* d, const TYPE_REAL* e, const int* m, \
                                 TYPE_REAL* w, const int* iblock, const int* isplit, const TYPE_REAL* orfac, \
                                 TYPE* Z, const int* iZ, const int* jZ, const int* descZ, \
                                 TYPE* work, const int* lwork, int* iwork, const int* liwork, \
                                 int* ifail, int* iclustr, TYPE_REAL* gap, int* info);

SCALAPACK_PSTEIN_DECL(psstein, PSSTEIN, float, float);
SCALAPACK_PSTEIN_DECL(pdstein, PDSTEIN, double, double);
SCALAPACK_PSTEIN_DECL(pcstein, PCSTEIN, lapack_complex_float, float);
SCALAPACK_PSTEIN_DECL(pzstein, PZSTEIN, lapack_complex_double, double);

#undef SCALAPACK_PSTEIN_DECL


#define SCALAPACK_PSTEDC_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* compz, const int* n, TYPE* d, TYPE* e, \
                                 TYPE* Q, const int* iq, const int* jq, const int* descQ, \
                                 TYPE* work, const int* lwork, int* iwork, const int* liwork, int* info);

SCALAPACK_PSTEDC_DECL(psstedc, PSSTEDC, float);
SCALAPACK_PSTEDC_DECL(pdstedc, PDSTEDC, double);

#undef SCALAPACK_PSTEDC_DECL

#ifdef __cplusplus
}
#endif

#endif // ROKKO_SCALAPACK_SCALAPACK_INTERFACE_H
