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

#ifndef ROKKO_XBLAS_LEVEL3_HPP
#define ROKKO_XBLAS_LEVEL3_HPP

#include <complex>
#include <stdexcept>
#include <cblas.h>
#include "util.hpp"

namespace rokko {
namespace xblas {

namespace {

template<typename T> struct gemm_dispatch;

template<>
struct gemm_dispatch<float> {
  template<typename MATRIX>
  static void gemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans_a,
                   enum CBLAS_TRANSPOSE trans_b, int m, int n, int k,
                   float alpha, MATRIX const& a, MATRIX const& b,
                   float beta, MATRIX& c) {
    cblas_sgemm(order, trans_a, trans_b, m, n, k,
                alpha, storage(a), lda(a), storage(b), lda(b),
                beta, storage(c), lda(c));
  }
};
  
template<>
struct gemm_dispatch<double> {
  template<typename MATRIX>
  static void gemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans_a,
                   enum CBLAS_TRANSPOSE trans_b, int m, int n, int k,
                   double alpha, MATRIX const& a, MATRIX const& b,
                   double beta, MATRIX& c) {
    cblas_dgemm(order, trans_a, trans_b, m, n, k,
                alpha, storage(a), lda(a), storage(b), lda(b),
                beta, storage(c), lda(c));
  }
};
  
template<>
struct gemm_dispatch<std::complex<float> > {
  template<typename MATRIX>
  static void gemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans_a,
                   enum CBLAS_TRANSPOSE trans_b, int m, int n, int k,
                   std::complex<float> const& alpha, MATRIX const& a,
                   MATRIX const& b,
                   std::complex<float> const& beta, MATRIX& c) {
    cblas_cgemm(order, trans_a, trans_b, m, n, k,
                &alpha, storage(a), lda(a), storage(b), lda(b),
                &beta, storage(c), lda(c));
  }
};
  
template<>
struct gemm_dispatch<std::complex<double> > {
  template<typename MATRIX>
  static void gemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE trans_a,
                   enum CBLAS_TRANSPOSE trans_b, int m, int n, int k,
                   std::complex<double> const& alpha, MATRIX const& a,
                   MATRIX const& b,
                   std::complex<double> const& beta, MATRIX& c) {
    cblas_zgemm(order, trans_a, trans_b, m, n, k,
                &alpha, storage(a), lda(a), storage(b), lda(b),
                &beta, storage(c), lda(c));
  }
};
  
}
  
template<typename MATRIX, typename T>
void gemm(enum CBLAS_TRANSPOSE trans_a, enum CBLAS_TRANSPOSE trans_b,
          T alpha, MATRIX const& a, MATRIX const& b, T beta, MATRIX& c) {
  if (util::op_cols(trans_a, a) != util::op_rows(trans_b, b))
    throw std::invalid_argument("matrix size mismatch");
  int m = util::op_rows(trans_a, a);
  int n = util::op_cols(trans_b, b);
  int k = util::op_cols(trans_a, a);
  gemm_dispatch<typename value_t<MATRIX>::type>::
    gemm((is_col_major(a) ? CblasColMajor : CblasRowMajor), trans_a, trans_b,
         m, n, k, alpha, a, b, beta, c);
}

} // namespace xblas
} // namespace rokko

#endif // ROKKO_XBLAS_LEVEL3_HPP
