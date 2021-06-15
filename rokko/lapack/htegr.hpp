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

template<typename T> struct htegr_dispatch;

template<>
struct htegr_dispatch<float> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htegr(int matrix_layout, char jobz, char range, lapack_int n, VECTOR& d, VECTOR& e,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
    return LAPACKE_sstegr(matrix_layout, jobz, range, n, storage(d), storage(e), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(isuppz));
  }
};

template<>
struct htegr_dispatch<double> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htegr(int matrix_layout, char jobz, char range, lapack_int n, VECTOR& d, VECTOR& e,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
    return LAPACKE_dstegr(matrix_layout, jobz, range, n, storage(d), storage(e), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(isuppz));
  }
};

template<>
struct htegr_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htegr(int matrix_layout, char jobz, char range, lapack_int n, VECTOR& d, VECTOR& e,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
    return LAPACKE_sstegr(matrix_layout, jobz, range, n, storage(d), storage(e), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(isuppz));
  }
};

template<>
struct htegr_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htegr(int matrix_layout, char jobz, char range, lapack_int n, VECTOR& d, VECTOR& e,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
    return LAPACKE_dstegr(matrix_layout, jobz, range, n, storage(d), storage(e), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(isuppz));
  }
};

} // end of anonymous namespace

template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(char jobz, char range, VECTOR& d, VECTOR& e,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  static_assert(std::is_same_v<real_t<MATRIX>, T>);
  static_assert(std::is_same_v<value_t<VECTOR>, T>);

  lapack_int n = size(d);
  if (size(e) != n)  // e(n) is used as internal workspace
    throw std::invalid_argument("vector e size mismatch");
  if (size(w) != n)
    throw std::invalid_argument("vector w size mismatch");
  return htegr_dispatch<value_t<MATRIX>>::
    htegr((is_col_major(z) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, isuppz);
}

// only eigenvalues (without jobz)
template<typename T, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(char range, VECTOR& d, VECTOR& e,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  static_assert(std::is_same_v<value_t<VECTOR>, T>);

  lapack_int n = size(d);
  constexpr null_matrix<T> z_null;
  return htegr_dispatch<T>::
    htegr(LAPACK_COL_MAJOR, 'N', range, n, d, e, vl, vu, il, iu, abstol, m, w, z_null, isuppz);
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(char range, VECTOR& d, VECTOR& e,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  return htegr('V', range, d, e, vl, vu, il, iu, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  constexpr real_t<VECTOR> real_zero = 0;
  return htegr('A', d, e, real_zero, real_zero, 0, 0, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  return htegr('V', d, e, vl, vu, 0, 0, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  constexpr real_t<VECTOR> real_zero = 0;
  return htegr('I', d, e, real_zero, real_zero, il, iu, abstol, m, w, z, isuppz);
}

// only eigenvalues, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  static_assert(std::is_same_v<value_t<VECTOR>, T>);
  constexpr real_t<VECTOR> real_zero = 0;

  return htegr('A', d, e, real_zero, real_zero, 0, 0, abstol, m, w, isuppz);
}

// only eigenvalues, use vl, vu (without jobz, range, il, iu)
template<typename T, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  static_assert(std::is_same_v<value_t<VECTOR>, T>);

  return htegr('V', d, e, vl, vu, 0, 0, abstol, m, w, isuppz);
}

// only eigenvalues, use il, iu (without jobz, range, vl, vu)
template<typename T, typename VECTOR, typename VECTOR_INT>
lapack_int htegr(VECTOR& d, VECTOR& e,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  constexpr real_t<VECTOR> real_zero = 0;
  return htegr('I', d, e, real_zero, real_zero, il, iu, abstol, m, w, isuppz);
}

ALIAS_TEMPLATE_FUNCTION(stegr, htegr);

} // end namespace lapack
} // end namespace rokko
