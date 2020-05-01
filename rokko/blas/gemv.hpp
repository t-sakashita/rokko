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

#ifndef ROKKO_BLAS_GEMV_HPP
#define ROKKO_BLAS_GEMV_HPP

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

#define BLAS_GEMV_IMPL(NAMES, TYPE) \
inline void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans, int m, int n, TYPE alpha, const TYPE * a, int lda, const TYPE * x, int inc_x, TYPE beta, TYPE * y, int inc_y) { \
  cblas_ ## NAMES (order, trans, m, n, alpha, a, lda, x, inc_x, beta, y, inc_y); \
}

BLAS_GEMV_IMPL(sgemv, float);
BLAS_GEMV_IMPL(dgemv, double);
  
#undef BLAS_GEMV_IMPL  
  
#define BLAS_GEMV_IMPL(NAMES, TYPE) \
inline void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans, int m, int n, TYPE alpha, const TYPE * a, int lda, const TYPE * x, int inc_x, TYPE beta, TYPE * y, int inc_y) { \
  cblas_ ## NAMES (order, trans, m, n, &alpha, lapack::complex_cast(a), lda, lapack::complex_cast(x), inc_x, &beta, lapack::complex_cast(y), inc_y); \
}

BLAS_GEMV_IMPL(cgemv, std::complex<float>);
BLAS_GEMV_IMPL(zgemv, std::complex<double>);
  
#undef BLAS_GEMV_IMPL  
  
template<typename MATRIX, typename VECTOR, typename T>
void gemv(CBLAS_TRANSPOSE trans, T alpha, const MATRIX& a, const VECTOR& x, int inc_x,
          T beta, VECTOR& y, int inc_y) {
  if (util::op_cols(trans, a) != size(x) / inc_x ||
      util::op_cols(trans, a) != size(y) / inc_y)
    throw std::invalid_argument("matrix/vector size mismatch");
  int m = util::op_rows(trans, a);
  int n = util::op_cols(trans, a);
  gemv((is_col_major(a) ? CblasColMajor : CblasRowMajor), trans, m, n, alpha, storage(a),
       ld(a), storage(x), inc_x, beta, storage(y), inc_y);
}

} // namespace blas
} // namespace rokko

#endif // ROKKO_BLAS_GEMV_HPP
