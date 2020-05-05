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

#ifndef ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP
#define ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_dsyevx.hpp>

namespace rokko {
namespace lapack {

// bisection (dsyevx) only eigenvalues
template<typename T, int MATRIX_MAJOR>
parameters diagonalize_bisection(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, T* eigvals,
                 rokko::parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<T>("abstol") < 0) {
      std::stringstream msg;
      msg << "lapack::diagonalize_bisection() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }

  parameters params_out = diagonalize_dsyevx(mat, eigvals, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (bisection)" << std::endl;
  }

  return params_out;
}


// bisection (dsyevx) eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR>
parameters diagonalize_bisection(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, T* eigvals,
				 Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
                 rokko::parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<T>("abstol") < 0) {
      std::stringstream msg;
      msg << "lapack::diagonalize_bisection() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }

  parameters params_out = diagonalize_dsyevx(mat, eigvals, eigvecs, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished dsyevx (bisection)" << std::endl;
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP
