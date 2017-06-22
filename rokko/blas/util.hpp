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

#ifndef ROKKO_BLAS_UTIL_HPP
#define ROKKO_BLAS_UTIL_HPP

#include <cblas.h>

namespace rokko {
namespace blas {
namespace util {

template<typename MATRIX>
int op_rows(enum CBLAS_TRANSPOSE trans, MATRIX const& a) {
  return (trans == CblasNoTrans) ? rows(a) : cols(a);
}
  
template<typename MATRIX>
int op_cols(enum CBLAS_TRANSPOSE trans, MATRIX const& a) {
  return (trans == CblasNoTrans) ? cols(a) : rows(a);
}
  
} // namespace util
} // namespace blas
} // namespace rokko

#endif // ROKKO_BLAS_UTIL_HPP
