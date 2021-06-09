/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <complex>
#ifdef I
# undef I
#endif
#include <stdexcept>
#include <cblas.h>
#include <rokko/blas/util.hpp>
#include <rokko/lapack/complex_cast.hpp>

namespace rokko {
namespace blas {

#define BLAS_GEMM_IMPL(NAMES, TYPE) \
  inline void gemm(CBLAS_ORDER order, CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b, int m, int n, int k, TYPE alpha, const TYPE * a, int lda_a, const TYPE * b, int lda_b, TYPE beta, TYPE * c, int lda_c) { \
  cblas_ ## NAMES (order, trans_a, trans_b, m, n, k, alpha, a, lda_a, b, lda_b, beta, c, lda_c); \
}
  
BLAS_GEMM_IMPL(sgemm, float);
BLAS_GEMM_IMPL(dgemm, double);
  
#undef BLAS_GEMM_IMPL  

#define BLAS_GEMM_IMPL(NAMES, TYPE) \
inline void gemm(CBLAS_ORDER order, CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b, int m, int n, int k, TYPE alpha, const TYPE * a, int lda_a, const TYPE * b, int lda_b, TYPE beta, TYPE * c, int lda_c) { \
  cblas_ ## NAMES (order, trans_a, trans_b, m, n, k, &alpha, lapack::complex_cast(a), lda_a, lapack::complex_cast(b), lda_b, &beta, lapack::complex_cast(c), lda_c); \
}
  
BLAS_GEMM_IMPL(cgemm, std::complex<float>);
BLAS_GEMM_IMPL(zgemm, std::complex<double>);
  
#undef BLAS_GEMM_IMPL  

template<typename MATRIX, typename T>
void gemm(CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b,
          T alpha, MATRIX const& a, MATRIX const& b, T beta, MATRIX& c) {
  if (util::op_cols(trans_a, a) != util::op_rows(trans_b, b))
    throw std::invalid_argument("matrix size mismatch");
  int m = util::op_rows(trans_a, a);
  int n = util::op_cols(trans_b, b);
  int k = util::op_cols(trans_a, a);
  gemm((is_col_major(a) ? CblasColMajor : CblasRowMajor), trans_a, trans_b,
       m, n, k, alpha, storage(a), ld(a), storage(b), ld(b), beta, storage(c), ld(c));
}

} // namespace blas
} // namespace rokko
