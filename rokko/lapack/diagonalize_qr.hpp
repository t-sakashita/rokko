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

#ifndef ROKKO_LAPACK_DIAGONALIZE_QR_HPP
#define ROKKO_LAPACK_DIAGONALIZE_QR_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_dsyevx.hpp>

namespace rokko {
namespace lapack {

// qr (dsyevx) only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_qr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
                 rokko::parameters params) {
  if (params.defined("abstol")) {
    double abstol = params.get<double>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  parameters params_out = diagonalize_dsyevx(mat, eigvals, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (qr)" << std::endl;
  }

  return params_out;
}


// qr (dsyevx) eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_qr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
                 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
                 rokko::parameters params) {
  if (params.defined("abstol")) {
    double abstol = params.get<double>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  parameters params_out = diagonalize_dsyevx(mat, eigvals, eigvecs, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (qr)" << std::endl;
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_QR_HPP