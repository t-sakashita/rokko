#ifndef ROKKO_PBLAS_H
#define ROKKO_PBLAS_H

#include <mpi.h>

extern "C" {
  void descinit_(int* desc, const int& m, const int& n, const int& mb, const int& nb,
                 const int& irsrc, const int& icsrc, const int& ixtxt, const int& lld, int& info);

  void pdgemm_(const char* TRANSA, const char* TRANSB,
	       const int& M, const int& N, const int& K,
	       const double& ALPHA,
	       double* A, const int& IA, const int& JA, int* DESCA,
	       double* B, const int& IB, const int& JB, int* DESCB,
	       const double& BETA,
	       double* C, const int& IC, const int& JC, int* DESCC);
}

namespace rokko {
namespace pblas {

template <class MATRIX>
void create_desc(const MATRIX& mat, int* desc)
{
  const int ZERO=0;
  int info;

  int lld = mat.m_local;
  //cout << "lld=" << lld << endl;
  if (lld == 0) lld = 1;
  descinit_(desc, mat.m_global, mat.n_global, mat.mb, mat.nb, ZERO, ZERO, mat.g.ictxt, lld, info);
}



template <class MATRIX>
void product(const MATRIX& matA, bool transA, const MATRIX& matB, bool transB, double alpha, double beta, MATRIX& matC)
{
  const int ZERO=0, ONE=1;
  int descA[9], descB[9], descC[9];
  int info;

  create_desc(matA, descA);
  create_desc(matB, descB);
  create_desc(matC, descC);

  char char_transA, char_transB;
  if (transA) char_transA = 'T';
  else char_transA = 'N';
  if (transB) char_transB = 'T';
  else char_transB = 'N';

  pdgemm_(&char_transA, &char_transB,
	  matA.m_global, matB.n_global, matA.n_global,
	  alpha,
	  matA.array, ONE, ONE, descA,
          matB.array, ONE, ONE, descB,
	  beta,
          matC.array, ONE, ONE, descC);
}

} // namespace pblas
} // namespace rokko

#endif // ROKKO_PBLAS_H
