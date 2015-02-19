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

#ifndef ROKKO_PBLAS_WRAP_H
#define ROKKO_PBLAS_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

void ROKKO_pdcopy(int N, const double* X, int IX, int JX, int* DESCX, int INCX,
                  double* Y, int IY, int JY, int *DESCY, int INCY);

double ROKKO_pddot(int N, const double* X, int IX, int JX, int* DESCX, int INCX,
                   const double* Y, int IY, int JY, int* DESCY, int INCY);

void ROKKO_pdgemv(char TRANS, int M, int N, double ALPHA,
                  const double* A, int IA, int JA, int* DESCA,
                  const double* X, int IX, int JX, int* DESCX, int INCX, double BETA,
                  double* Y, int IY, int JY, int* DESCY, int INCY);

void ROKKO_pdgemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
                  const double* A, int IA, int JA, int* DESCA,
                  const double* B, int IB, int JB, int* DESCB,
                  double BETA, double* C, int IC, int JC, int* DESCC);
  
#ifdef __cplusplus
}
#endif

#endif // ROKKO_PBLAS_WRAP_H
