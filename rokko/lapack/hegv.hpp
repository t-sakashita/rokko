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

#pragma once

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

template<typename T> struct hegv_dispatch;

template<>
struct hegv_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegv(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                         MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_ssygv(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegv_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegv(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                         MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_dsygv(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegv_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegv(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                         MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_chegv(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

template<>
struct hegv_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int hegv(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                         MATRIX& a, MATRIX& b, VECTOR& w) {
    return LAPACKE_zhegv(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), storage(w));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
lapack_int hegv(lapack_int itype, char jobz, char uplo, MATRIX& a, MATRIX& b, VECTOR& w) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR>>);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return hegv_dispatch<value_t<MATRIX>>::
    hegv((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), itype, jobz, uplo, n, a, b, w);
}

ALIAS_TEMPLATE_FUNCTION(sygv, hegv);

} // end namespace lapack
} // end namespace rokko
