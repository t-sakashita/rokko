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

#ifndef ROKKO_LAPACK_HTEBZ_HPP
#define ROKKO_LAPACK_HTEBZ_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct htebz_dispatch;

template<>
struct htebz_dispatch<float> {
  template<typename VECTOR, typename VECTOR_INT>
  static lapack_int htebz(char range, char order, lapack_int n, float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
    return LAPACKE_sstebz(range, order, n, vl, vu, il, iu, abstol, storage(d), storage(e), &m, &nsplit, storage(w), storage(iblock), storage(isplit));
  }
};

template<>
struct htebz_dispatch<double> {
  template<typename VECTOR, typename VECTOR_INT>
  static lapack_int htebz(char range, char order, lapack_int n, double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
    return LAPACKE_dstebz(range, order, n, vl, vu, il, iu, abstol, storage(d), storage(e), &m, &nsplit, storage(w), storage(iblock), storage(isplit));
  }
};

template<>
struct htebz_dispatch<std::complex<float>> {
  template<typename VECTOR, typename VECTOR_INT>
  static lapack_int htebz(char range, char order, lapack_int n, float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
    return LAPACKE_cstebz(range, order, n, vl, vu, il, iu, abstol, storage(d), storage(e), &m, &nsplit, storage(w), storage(iblock), storage(isplit));
  }
};

template<>
struct htebz_dispatch<std::complex<double>> {
  template<typename VECTOR, typename VECTOR_INT>
  static lapack_int htebz(char range, char order, lapack_int n, double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
    return LAPACKE_zstebz(range, order, n, vl, vu, il, iu, abstol, storage(d), storage(e), &m, &nsplit, storage(w), storage(iblock), storage(isplit));
  }
};

} // end of anonymous namespace

template<typename VECTOR, typename VECTOR_INT>
lapack_int htebz(char range, char order, double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                 VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
  lapack_int n = size(d);
  if (size(e) != (n-1))
    throw std::invalid_argument("vector e size mismatch");
  if (size(w) != n)
    throw std::invalid_argument("vector w size mismatch");
  return htebz_dispatch<value_t<VECTOR>>::
    htebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit);
}

template<typename VECTOR, typename VECTOR_INT>
lapack_int stebz(char range, char order, double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                 VECTOR& d, VECTOR& e, lapack_int& m, lapack_int& nsplit, VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
  return htebz(range, order, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit);
};

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HTEBZ_HPP
