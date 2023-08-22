/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_syevx.hpp>

namespace rokko {
namespace lapack {

// qr (dsyevx) only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_qr(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
                 rokko::parameters params) {
  if (params.defined("abstol")) {
    const auto abstol = params.get<real_t<T>>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  const auto params_out = diagonalize_syevx(mat, eigvals, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (qr)" << std::endl;
  }

  return params_out;
}


// qr (dsyevx) eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_qr(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
                 Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
                 rokko::parameters params) {
  if (params.defined("abstol")) {
    const auto abstol = params.get<real_t<T>>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  const auto params_out = diagonalize_syevx(mat, eigvals, eigvecs, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (qr)" << std::endl;
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko
