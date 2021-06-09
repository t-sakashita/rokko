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

#pragma once

#include <rokko/cscalapack.h>

namespace rokko {
namespace scalapack {

namespace {

inline void pelset_dispatch(double* A, int ia, int ja, const int* descA, double alpha) {
  cscalapack_pdelset(A, ia, ja, descA, alpha);
}
  
inline void pelset_dispatch(float* A, int ia, int ja, const int* descA, float alpha) {
  cscalapack_pselset(A, ia, ja, descA, alpha);
}
  
}

template<typename MATRIX>
void pelset(MATRIX& A, int ia, int ja, typename MATRIX::value_type alpha) {
  pelset_dispatch(A.get_array_pointer(), ia, ja, A.get_blacs_descriptor(), alpha);
}

} // end namespace scalapack
} // end namespace rokko
