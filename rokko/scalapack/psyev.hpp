/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_PSYEV_HPP
#define ROKKO_SCALAPACK_PSYEV_HPP

#include <rokko/cscalapack.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

namespace {

inline int psyev_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                          double* w, double* Z, int iz, int jz, const int* descZ,
                          double* work, int lwork) {
  return cscalapack_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, work, lwork);
}
  
inline int psyev_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                          double* w, double* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pdsyev(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ);
}
  
}
  
template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, int n, MATRIX& a, const int* descA, VECTOR& w,
          MATRIX& z, const int* descZ) {
  return psyev_dispatch(jobz, uplo, n, storage(a), 0, 0, descA, storage(w),
                        storage(z), 0, 0, descZ);
}

template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, int n, MATRIX& a, const int* descA, VECTOR& w,
          MATRIX& z, const int* descZ, VECTOR& work) {
  return psyev_dispatch(jobz, uplo, n, storage(a), 0, 0, descA, storage(w),
                        storage(z), 0, 0, descZ, storage(work), work.size());
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEV_HPP
