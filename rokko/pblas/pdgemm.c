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

#include <rokko/pblas/pblas.h>

void ROKKO_pdgemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
                  const double* A, int IA, int JA, int* DESCA,
                  const double* B, int IB, int JB, int* DESCB,
                  double BETA, double* C, int IC, int JC, int* DESCC) {
  pdgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA, B, &IB, &JB, DESCB,
          &BETA, C, &IC, &JC, DESCC);
}
