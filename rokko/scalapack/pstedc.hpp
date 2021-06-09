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

#ifndef ROKKO_SCALAPACK_PSTEDC_HPP
#define ROKKO_SCALAPACK_PSTEDC_HPP

#include <lapacke.h>
#include <rokko/cscalapack.h>
#include <rokko/eigen3.hpp>
#include <rokko/lapack/storage.hpp>
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

namespace {

inline int pstedc_dispatch(char compz, int n, double* d, double* e,
                           double* Q, int iq, int jq, const int* descQ) {
  return cscalapack_pdstedc(compz, n, d, e,
                            Q, iq, jq, descQ);
}

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
int pstedc(char compz, VECTOR& d, VECTOR& e, MATRIX& z) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR>>);

  lapack_int n = size(d);
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return pstedc_dispatch(compz, n, storage(d), storage(e),
                         z.get_array_pointer(), 0, 0, descZ);
}

// eigenvalues & eigenvectors (without compz)
template<typename MATRIX, typename VECTOR>
int pstedc(MATRIX& a, VECTOR& w, MATRIX& z) {
  return pstedc('I', a, w, z);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSTEDC_HPP
