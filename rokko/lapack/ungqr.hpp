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

#ifndef ROKKO_LAPACK_UNGQR_HPP
#define ROKKO_LAPACK_UNGQR_HPP

#include <complex>
#include <stdexcept>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <lapacke.h>
#undef I
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct ungqr_dispatch;
  
template<>
struct ungqr_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_sorgqr(matrix_layout, m, n, k, storage(a), lda(a), storage(tau));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_sorgqr(matrix_layout, m, n, k, storage(a), lda(a), storage(tau),
                          storage(work), size(work));
  }
};

template<>
struct ungqr_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_dorgqr(matrix_layout, m, n, k, storage(a), lda(a), storage(tau));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_dorgqr(matrix_layout, m, n, k, storage(a), lda(a), storage(tau),
                          storage(work), size(work));
  }
};

template<>
struct ungqr_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_cungqr(matrix_layout, m, n, k, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_cungqr(matrix_layout, m, n, k, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)), complex_cast(storage(work)),
                          size(work));
  }
};

template<>
struct ungqr_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_zungqr(matrix_layout, m, n, k, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int ungqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_zungqr(matrix_layout, m, n, k, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)), complex_cast(storage(work)),
                          size(work));
  }
};

}
  
template<typename MATRIX, typename VECTOR>
lapack_int ungqr(lapack_int k, MATRIX& a, VECTOR const& tau) {
  BOOST_STATIC_ASSERT(std::is_same<value_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  if (size(tau) != k)
    throw std::invalid_argument("vector tau size mismatch");
  return ungqr_dispatch<value_t<MATRIX>>
    ::ungqr((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), m, n, k, a, tau);
}

template<typename MATRIX, typename VECTOR>
lapack_int ungqr(lapack_int k, MATRIX& a, VECTOR const& tau, VECTOR& work) {
  BOOST_STATIC_ASSERT(std::is_same<value_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  if (size(tau) != k)
    throw std::invalid_argument("vector tau size mismatch");
  if (size(work) < std::max(1, n))
    throw std::invalid_argument("vector work size mismatch");
  return ungqr_dispatch<value_t<MATRIX>>
    ::ungqr((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), m, n, k, a, tau, work);
}

template<typename MATRIX, typename VECTOR>
lapack_int orgqr(lapack_int k, MATRIX& a, VECTOR const& tau) {
  return ungqr(k, a, tau);
}
  
template<typename MATRIX, typename VECTOR>
lapack_int orgqr(lapack_int k, MATRIX& a, VECTOR const& tau, VECTOR& work) {
  return ungqr(k, a, tau, work);
};
  
} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_UNGQR_HPP
