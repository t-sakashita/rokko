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

#ifndef ROKKO_LAPACK_GETRF_HPP
#define ROKKO_LAPACK_GETRF_HPP

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

template<typename T> struct getrf_dispatch;
  
template<>
struct getrf_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_sgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

template<>
struct getrf_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_dgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

template<>
struct getrf_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_cgetrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          storage(ipiv));
  }
};

template<>
struct getrf_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_zgetrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          storage(ipiv));
  }
};

}

template<typename MATRIX, typename VECTOR>
lapack_int getrf(MATRIX& a, VECTOR& ipiv) {
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  BOOST_STATIC_ASSERT(std::is_same<value_t<VECTOR>, lapack_int>::value);
  if (size(ipiv) < std::min(m, n))
    throw std::invalid_argument("vector ipiv size mismatch");
  return getrf_dispatch<value_t<MATRIX>>
    ::getrf((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
            m, n, a, ipiv);
}

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_GETRF_HPP
