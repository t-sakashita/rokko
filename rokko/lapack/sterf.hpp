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

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct sterf_dispatch;

template<>
struct sterf_dispatch<float> {
  template<typename VECTOR>
  static lapack_int sterf(lapack_int n,
                          VECTOR& d, VECTOR& e) {
    return LAPACKE_ssterf(n, storage(d), storage(e));
  }
};

template<>
struct sterf_dispatch<double> {
  template<typename VECTOR>
  static lapack_int sterf(lapack_int n,
                          VECTOR& d, VECTOR& e) {
    return LAPACKE_dsterf(n, storage(d), storage(e));
  }
};

} // end of anonymous namespace

template<typename VECTOR>
lapack_int sterf(VECTOR& d, VECTOR& e) {
  lapack_int n = size(d);
  if (size(e) != (n-1))
    throw std::invalid_argument("vector e size mismatch");
  return sterf_dispatch<value_t<VECTOR>>::
    sterf(n, d, e);
}

} // end namespace lapack
} // end namespace rokko
