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

#ifndef ROKKO_SCALAPACK_PLANGE_HPP
#define ROKKO_SCALAPACK_PLANGE_HPP

#include <rokko/cscalapack.h>

namespace rokko {
namespace scalapack {

namespace {

inline double plange_dispatch(char norm, int m, int n, const double* A, const int* descA) {
  return cscalapack_pdlange(norm, m, n, A, descA);
}

inline float plange_dispatch(char norm, int m, int n, const float* A, const int* descA) {
  return cscalapack_pslange(norm, m, n, A, descA);
}
  
}

template<typename MATRIX>
typename MATRIX::value_type plange(char norm, const MATRIX& A) {
  return plange_dispatch(norm, A.get_m_global(), A.get_n_global(), A.get_array_pointer(),
                         A.get_mapping().get_blacs_descriptor());
}

template<typename MATRIX, typename VECTOR>
typename MATRIX::value_type plange(char norm, const MATRIX& A, VECTOR& work) {
  return plange_dispatch(norm, A.get_m_global(), A.get_n_global(), A.get_array_pointer(),
                         A.get_mapping().get_blacs_descriptor(), work);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PLANGE_HPP
