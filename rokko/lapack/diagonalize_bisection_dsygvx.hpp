/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
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

namespace rokko {
namespace lapack {

// only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_bisection(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			  VEC& eigvals,
			  rokko::parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<real_t<T>>("abstol") < 0) {
      std::stringstream msg;
      msg << "lapack::diagonalize_bisection_sygvx() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use sygvx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }

  return diagonalize_sygvx(mata, matb, eigvals, params);
}


// eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_bisection(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			  VEC& eigvals,
			  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			  parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<real_t<T>>("abstol") < 0) {
      std::stringstream msg;
      msg << "lapack::diagonalize_bisection_sygvx() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use sygvx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }

  return diagonalize_sygvx(mata, matb, eigvals, eigvecs, params);
}


} // namespace lapack
} // namespace rokko
