/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_PSTEBZ_HPP
#define ROKKO_SCALAPACK_PSTEBZ_HPP

#include <lapacke.h>
#include <rokko/cscalapack.h>
#include <rokko/eigen3.hpp>
#include <rokko/lapack/storage.hpp>
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

using rokko::lapack::storage;

namespace {

inline int pstebz_dispatch(int ictxt, char range, char order, int n,
                           double vl, double vu, int il, int iu,
                           double abstol, const double* d, const double* e, int& m, int& nsplit,
                           double* w, int* iblock, int* isplit) {
  return cscalapack_pdstebz(ictxt, range, order, n,
                            vl, vu, il, iu,
                            abstol, d, e, &m, &nsplit,
                            w, iblock, isplit);
}

} // end of anonymous namespace

template<typename T, typename VECTOR, typename VECTOR_INT>
int pstebz(int ictxt, char range, char order,
           T vl, T vu, int il, int iu,
           T abstol, const VECTOR& d, const VECTOR& e, int& m, int& nsplit,
           VECTOR& w, VECTOR_INT& iblock, VECTOR_INT& isplit) {
  static_assert(std::is_same<value_t<VECTOR>, T>::value, "");
  lapack_int n = size(d);

  return pstebz_dispatch(ictxt, range, order, n,
                         vl, vu, il, iu,
                         abstol, storage(d), storage(e), m, nsplit,
                         storage(w), storage(iblock), storage(isplit));
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSTEBZ_HPP
