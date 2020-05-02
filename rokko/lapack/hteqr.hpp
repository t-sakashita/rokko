/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_HTEQR_HPP
#define ROKKO_LAPACK_HTEQR_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct hteqr_dispatch;

template<>
struct hteqr_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hteqr(int matrix_layout, char compz, lapack_int n,
                          VECTOR& d, VECTOR& e, MATRIX& z) {
    return LAPACKE_ssteqr(matrix_layout, compz, n, storage(d), storage(e), storage(z), ld(z));
  }
};

template<>
struct hteqr_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hteqr(int matrix_layout, char compz, lapack_int n,
                          VECTOR& d, VECTOR& e, MATRIX& z) {
    return LAPACKE_dsteqr(matrix_layout, compz, n, storage(d), storage(e), storage(z), ld(z));
  }
};

template<>
struct hteqr_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hteqr(int matrix_layout, char compz, lapack_int n,
                          VECTOR& d, VECTOR& e, MATRIX& z) {
    return LAPACKE_csteqr(matrix_layout, compz, n, storage(d), storage(e), complex_cast(storage(z)), ld(z));
  }
};

template<>
struct hteqr_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hteqr(int matrix_layout, char compz, lapack_int n,
                          VECTOR& d, VECTOR& e, MATRIX& z) {
    return LAPACKE_zsteqr(matrix_layout, compz, n, storage(d), storage(e), complex_cast(storage(z)), ld(z));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
lapack_int hteqr(char compz, VECTOR& d, VECTOR& e, MATRIX& z) {
  BOOST_STATIC_ASSERT(std::is_same<norm_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int n = rows(z);
  if (rows(z) != cols(z))
    throw std::invalid_argument("matrix Z size mismatch");
  if (size(d) != n)
    throw std::invalid_argument("vector d size mismatch");
  if (size(e) != (n-1))
    throw std::invalid_argument("vector e size mismatch");
  return hteqr_dispatch<value_t<MATRIX>>::
    hteqr((is_col_major(z) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), compz, n, d, e, z);
}

ALIAS_TEMPLATE_FUNCTION(steqr, hteqr);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HTEQR_HPP
