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

#ifndef ROKKO_LAPACK_HEEVD_HPP
#define ROKKO_LAPACK_HEEVD_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/lapack/storage.hpp>
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct heevd_dispatch;

template<>
struct heevd_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heevd(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_ssyevd(matrix_layout, jobz, uplo, n, storage(a), ld(a), storage(w));
  }
};

template<>
struct heevd_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heevd(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, storage(a), ld(a), storage(w));
  }
};

template<>
struct heevd_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heevd(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_cheevd(matrix_layout, jobz, uplo, n, complex_cast(storage(a)), ld(a),
                         storage(w));
  }
};

template<>
struct heevd_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heevd(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_zheevd(matrix_layout, jobz, uplo, n, complex_cast(storage(a)), ld(a),
                         storage(w));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
lapack_int heevd(char jobz, char uplo, MATRIX& a, VECTOR& w) {
  BOOST_STATIC_ASSERT(std::is_same<norm_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return heevd_dispatch<value_t<MATRIX>>::
    heevd((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, uplo, n, a, w);
}

ALIAS_TEMPLATE_FUNCTION(syevd, heevd);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HEEVD_HPP
