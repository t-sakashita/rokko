#ifndef ROKKO_PBLAS_H
#define ROKKO_PBLAS_H

#include <rokko/scalapack/blacs.hpp>

extern "C" {
  void descinit_(int* desc, const int& m, const int& n, const int& mb, const int& nb,
                 const int& irsrc, const int& icsrc, const int& ixtxt, const int& lld, int& info);

  void pdgemm_(const char* TRANSA, const char* TRANSB,
	       const int& M, const int& N, const int& K,
	       const double& ALPHA,
	       const double* A, const int& IA, const int& JA, int* DESCA,
	       const double* B, const int& IB, const int& JB, int* DESCB,
	       const double& BETA,
	       double* C, const int& IC, const int& JC, int* DESCC);
}

namespace rokko {
namespace pblas {

template <class MATRIX>
void create_desc(const MATRIX& mat, const int ictxt, int* desc) {
  const int ZERO=0;
  int info;
  int lld = std::max(mat.get_m_local(), 1);
  descinit_(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(), mat.get_nb(), ZERO, ZERO,
            ictxt, lld, info);
}

template <class MATRIX>
void product(const MATRIX& matA, bool transA, const MATRIX& matB, bool transB, double alpha,
             double beta, MATRIX& matC) {
  const int MINUS_ONE = -1, ZERO=0, ONE=1;
  int ictxt;
  blacs_get_(MINUS_ONE, ZERO, ictxt);

  char char_grid_major = (matA.get_grid().is_row_major() ? 'R' : 'C');
  blacs_gridinit_(ictxt, &char_grid_major, matA.get_grid().get_nprow(), matA.get_grid().get_npcol());

  int descA[9], descB[9], descC[9];
  create_desc(matA, ictxt, descA);
  create_desc(matB, ictxt, descB);
  create_desc(matC, ictxt, descC);

  char char_transA = (transA ? 'T' : 'N');
  char char_transB = (transB ? 'T' : 'N');
  pdgemm_(&char_transA, &char_transB,
	  matA.get_m_global(), matB.get_n_global(), matA.get_n_global(),
	  alpha,
	  matA.get_array_pointer(), ONE, ONE, descA,
          matB.get_array_pointer(), ONE, ONE, descB,
	  beta,
          matC.get_array_pointer(), ONE, ONE, descC);

  blacs_gridexit_(&ictxt);
}

} // namespace pblas
} // namespace rokko

#endif // ROKKO_PBLAS_H
