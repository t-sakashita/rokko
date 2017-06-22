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

#ifndef ROKKO_LAPACK_COMPLEX_CAST_HPP
#define ROKKO_LAPACK_COMPLEX_CAST_HPP

#include <complex>

namespace rokko {
namespace lapack {

lapack_complex_float* complex_cast(float* data) {
  return reinterpret_cast<lapack_complex_float*>(data);
}
  
const lapack_complex_float* complex_cast(const float* data) {
  return reinterpret_cast<const lapack_complex_float*>(data);
}
  
lapack_complex_double* complex_cast(double* data) {
  return reinterpret_cast<lapack_complex_double*>(data);
}
  
const lapack_complex_double* complex_cast(const double* data) {
  return reinterpret_cast<const lapack_complex_double*>(data);
}
  
} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_COMPLEX_CAST_HPP
