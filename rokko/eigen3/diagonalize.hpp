/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN3_DIAGONALIZE_H
#define ROKKO_EIGEN3_DIAGONALIZE_H

#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <Eigen/Dense>

namespace rokko {
namespace eigen3 {

template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int info = 0;
  timer.stop(timer_id::diagonalize_initialize);
  timer.start(timer_id::diagonalize_diagonalize);
  Eigen::SelfAdjointEigenSolver<typename localized_matrix<MATRIX_MAJOR>::super_type> ES(mat);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  eigvals = ES.eigenvalues();
  eigvecs = ES.eigenvectors();
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, std::vector<double>& eigvals_in,
		localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int info = 0;
  int dim = mat.rows();
  if (eigvals_in.size() < dim) eigvals_in.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  Eigen::Map<Eigen::VectorXd> eigvals(&eigvals_in[0], eigvals_in.size());
  timer.start(timer_id::diagonalize_diagonalize);
  Eigen::SelfAdjointEigenSolver<typename localized_matrix<MATRIX_MAJOR>::super_type> ES(mat);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  eigvals = ES.eigenvalues();
  eigvecs = ES.eigenvectors();
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

} // namespace eigen3
} // namespace rokko

#endif // ROKKO_EIGEN3_DIAGONALIZE_H
