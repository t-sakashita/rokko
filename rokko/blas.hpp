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

#ifndef ROKKO_BLAS_HPP
#define ROKKO_BLAS_HPP

#include <cblas.h>

namespace rokko {
namespace blas {

// level 2

template<typename MATRIX, typename VECTOR, typename T>
void gemv(enum CBLAS_TRANSPOSE trans,
          T alpha, MATRIX const& a, VECTOR const& x, int inc_x,
          T beta, VECTOR& y, int inc_y);

// level 3

template<typename MATRIX, typename T>
void gemm(enum CBLAS_TRANSPOSE trans_a, enum CBLAS_TRANSPOSE trans_b,
          T alpha, MATRIX const& a, MATRIX const& b, T beta, MATRIX& c);

} // namespace blas
} // namespace rokko

#include "blas/level2.hpp"
#include "blas/level3.hpp"

#endif // ROKKO_BLAS_HPP
