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

#ifndef ROKKO_SCALAPACK_WRAP_H
#define ROKKO_SCALAPACK_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

void ROKKO_pdelget(char scope, char top, double* alpha, const double* A, int ia, int ja,
                   const int* descA);

void ROKKO_pdelset(double* A, int ia, int ja, const int* descA, double alpha);

double ROKKO_pdlamch(int icnt, char cmch);

void ROKKO_pdlaprnt(int m, int n, const double* A, int ia, int ja, const int* descA, int irprnt,
                    int icprnt, const char* cmatnm, int nout, double* work);

int ROKKO_pdsyev_work(char jobz, char uplo, int n,
		      double* A, int ia, int ja, const int* descA,
		      double* w, double* Z, int iz, int jz, const int* descZ,
		      double* work, int lwork);

int ROKKO_pdsyevd_work(char jobz, char uplo, int n,
		       double* A, int ia, int ja, const int* descA,
		       double* w, double* Z, int iz, int jz, const int* descZ,
		       double* work, int lwork, int* iwork, int liwork);

int ROKKO_pdsyev(char jobz, char uplo, int n,
		 double* A, int ia, int ja, const int* descA,
		 double* w, double* Z, int iz, int jz, const int* descZ);

int ROKKO_pdsyevd(char jobz, char uplo, int n,
                  double* A, int ia, int ja, const int* descA,
                  double* w, double* Z, int iz, int jz, const int* descZ);

int ROKKO_pdsyevr(char jobz, char uplo, int n,
                  double* A, int ia, int ja, const int* descA,
                  double vl, double vu, int il, int iu,
                  int* m, int* nz, double* w,
                  double* Z, int iz, int jz, const int* descZ,
                  double* work, int lwork, int* iwork, int liwork);

int ROKKO_pdsyevx(char jobz, char range, char uplo, int n,
                  double* A, int iA, int jA, const int* descA,
                  double vl, double vu, int il, int iu,
                  double abstol, int* m, int* nZ, double* w, double orfac,
                  double* Z, int iZ, int jZ, const int* descZ,
                  double* work, int lwork, int* iwork, int liwork,
                  int* ifail, int* iclustr, double* gap);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_SCALAPACK_WRAP_H
