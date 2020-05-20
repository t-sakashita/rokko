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

#ifndef ROKKO_SCALAPACK_PSTEIN_HPP
#define ROKKO_SCALAPACK_PSTEIN_HPP

#include <lapacke.h>
#include <rokko/cscalapack.h>
#include <rokko/eigen3.hpp>
#include <rokko/lapack/storage.hpp>
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

using rokko::lapack::storage;

namespace {

inline int pstein_dispatch(int n, const double* d, const double* e, int m,
                           double* w, const int* iblock, const int* isplit, double orfac,
                           double* Z, const int* iZ, const int* jZ, const int* descZ,
                           int* ifail, int* iclustr, double* gap) {
  return cscalapack_pdstein(n, d, e, m,
                            w, iblock, isplit, orfac,
                            Z, iZ, jZ, descZ,
                            ifail, iclustr, gap);
}

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
int pstein(const VECTOR& d, const VECTOR& e, int& m,
           VECTOR& w, const VECTOR_INT& iblock, const VECTOR_INT& isplit, double orfac,
           MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR& gap) {
  static_assert(std::is_same<norm_t<MATRIX>, value_t<VECTOR>>::value, "");
  lapack_int n = size(d);
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();

  return pstein_dispatch(n, storage(d), storage(e), m,
                         storage(w), storage(iblock), storage(isplit), orfac,
                         z.get_array_pointer(), 0, 0, descZ,
                         storage(ifail), storage(iclustr), storage(gap));
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSTEIN_HPP
