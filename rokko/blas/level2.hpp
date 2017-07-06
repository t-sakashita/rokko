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

#ifndef ROKKO_BLAS_LEVEL2_HPP
#define ROKKO_BLAS_LEVEL2_HPP

#include <rokko/traits/value_t.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <complex>
#include <stdexcept>
#include <cblas.h>
#include "util.hpp"

namespace rokko {
namespace blas {

namespace {

template<typename T> struct gemv_dispatch;

template<>
struct gemv_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                   int m, int n, float alpha, MATRIX const& a,
                   VECTOR const& x, int inc_x, 
                   float beta, VECTOR& y, int inc_y) {
    cblas_dgemv(order, trans, m, n, alpha, storage(a), lda(a),
                storage(x), inc_x, beta, storage(y), inc_y);
  }
};
  
template<>
struct gemv_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                   int m, int n, double alpha, MATRIX const& a,
                   VECTOR const& x, int inc_x,
                   double beta, VECTOR& y, int inc_y) {
    cblas_dgemv(order, trans, m, n, alpha, storage(a), lda(a),
                storage(x), inc_x, beta, storage(y), inc_y);
  }
};
  
template<>
struct gemv_dispatch<std::complex<float> > {
  template<typename MATRIX, typename VECTOR>
  static void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                   int m, int n, std::complex<float> const& alpha, MATRIX const& a,
                   VECTOR const& x, int inc_x, 
                   std::complex<float> const& beta, VECTOR& y, int inc_y) {
    cblas_cgemv(order, trans, m, n, &alpha, storage(a), lda(a),
                storage(x), inc_x, &beta, storage(y), inc_y);
  }
};
  
template<>
struct gemv_dispatch<std::complex<double> > {
  template<typename MATRIX, typename VECTOR>
  static void gemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                   int m, int n, std::complex<double> const& alpha, MATRIX const& a,
                   VECTOR const& x, int inc_x, 
                   std::complex<double> const& beta, VECTOR& y, int inc_y) {
    cblas_zgemv(order, trans, m, n, &alpha, storage(a), lda(a),
                storage(x), inc_x, &beta, storage(y), inc_y);
  }
};
  
}
  
template<typename MATRIX, typename VECTOR, typename T>
void gemv(CBLAS_TRANSPOSE trans,
          T alpha, MATRIX const& a, VECTOR const& x, int inc_x,
          T beta, VECTOR& y, int inc_y) {
  BOOST_STATIC_ASSERT(boost::is_same<typename value_t<MATRIX>::type,
                      typename value_t<VECTOR>::type>::value);
  if (util::op_cols(trans, a) != size(x) / inc_x ||
      util::op_cols(trans, a) != size(y) / inc_y)
    throw std::invalid_argument("matrix/vector size mismatch");
  int m = util::op_rows(trans, a);
  int n = util::op_cols(trans, a);
  gemv_dispatch<typename value_t<MATRIX>::type>
    ::gemv((is_col_major(a) ? CblasColMajor : CblasRowMajor), trans, m, n, alpha, a,
           x, inc_x, beta, y, inc_y);
}

} // namespace blas
} // namespace rokko

#endif // ROKKO_BLAS_LEVEL2_HPP
