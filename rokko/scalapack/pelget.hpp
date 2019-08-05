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

#ifndef ROKKO_SCALAPACK_PELGET_HPP
#define ROKKO_SCALAPACK_PELGET_HPP

#include <rokko/cscalapack.h>

namespace rokko {
namespace scalapack {

namespace {

inline double pelget_dispatch(char scope, char top, const double* A, int ia, int ja,
                              const int* descA) {
  return cscalapack_pdelget(scope, top, A, ia, ja, descA);
}
  
inline float pelget_dispatch(char scope, char top, const float* A, int ia, int ja,
                             const int* descA) {
  return cscalapack_pselget(scope, top, A, ia, ja, descA);
}
  
}

template<typename MATRIX>
typename MATRIX::value_type pelget(char scope, char top, const MATRIX& A, int ia, int ja)  {
  return pelget_dispatch(scope, top, A.get_array_pointer(), ia, ja, A.get_blacs_descriptor());
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PELGET_HPP
