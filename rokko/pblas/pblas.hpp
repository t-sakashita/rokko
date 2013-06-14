#ifndef ROKKO_PBLAS_H
#define ROKKO_PBLAS_H

#include <rokko/scalapack/blacs.hpp>

extern "C" {

void pdgemm_(const char* TRANSA, const char* TRANSB,
             const int* M, const int* N, const int* K,
             const double* ALPHA,
             const double* A, const int* IA, const int* JA, int* DESCA,
             const double* B, const int* IB, const int* JB, int* DESCB,
             const double* BETA, double* C, const int* IC, const int* JC, int* DESCC);

inline void ROKKO_pdgemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
                         const double* A, int IA, int JA, int* DESCA,
                         const double* B, int IB, int JB, int* DESCB,
                         double BETA, double* C, int IC, int JC, int* DESCC) {
  pdgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA, B, &IB, &JB, DESCB,
          &BETA, C, &IC, &JC, DESCC);
}
  
}

#endif // ROKKO_PBLAS_H
