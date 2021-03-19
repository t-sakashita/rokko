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

#ifndef ROKKO_LAPACK_HEGVD_HPP
#define ROKKO_LAPACK_HEGVD_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/lapack/storage.hpp>
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct hegvd_dispatch;

template<>
struct hegvd_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegvd(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_ssygvd(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegvd_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegvd(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_dsygvd(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegvd_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegvd(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_chegvd(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegvd_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegvd(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_zhegvd(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
lapack_int hegvd(lapack_int itype, char jobz, char uplo, MATRIX& a, MATRIX& b, VECTOR& w) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return hegvd_dispatch<value_t<MATRIX>>::
    hegvd((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), itype, jobz, uplo, n, a, b, w);
}

ALIAS_TEMPLATE_FUNCTION(sygvd, hegvd);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HEGVD_HPP
