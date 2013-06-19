/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PBLAS_PBLAS_H
#define ROKKO_PBLAS_PBLAS_H

#include <rokko/scalapack/blacs.hpp>

#ifdef __cplusplus
extern "C" {
#endif

void pdgemm_(const char* TRANSA, const char* TRANSB,
             const int* M, const int* N, const int* K,
             const double* ALPHA,
             const double* A, const int* IA, const int* JA, int* DESCA,
             const double* B, const int* IB, const int* JB, int* DESCB,
             const double* BETA, double* C, const int* IC, const int* JC, int* DESCC);

void ROKKO_pdgemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
                  const double* A, int IA, int JA, int* DESCA,
                  const double* B, int IB, int JB, int* DESCB,
                  double BETA, double* C, int IC, int JC, int* DESCC);
  
#ifdef __cplusplus
}
#endif

#endif // ROKKO_PBLAS_PBLAS_H
