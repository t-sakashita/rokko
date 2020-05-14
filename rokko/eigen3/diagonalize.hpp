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

#ifndef ROKKO_EIGEN3_DIAGONALIZE_HPP
#define ROKKO_EIGEN3_DIAGONALIZE_HPP

#include <rokko/eigen3.hpp>
#include <rokko/traits/norm_t.hpp>
#include <rokko/parameters.hpp>

namespace rokko {
namespace eigen3 {

char get_matrix_part(rokko::parameters const& params) {
  std::string matrix_part;
  if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");
  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  if ((matrix_part[0] == 'u') || (matrix_part[0] == 'U'))
    return 'U';
  if ((matrix_part[0] == 'l') || (matrix_part[0] == 'L'))
    return 'L';
  else
    return '\0';
}

// only eigenvalues
template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, int SIZE>
parameters diagonalize(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, Eigen::Vector<norm_t<T>, SIZE>& eigvals,
		       rokko::parameters const& params) {
  parameters params_out;
  if (get_matrix_part(params) == 'U') {
    throw std::invalid_argument("eigen3::diagonalize() : Eigen3's SelfAdjointEigenSolver does not support upper part of matrix.");
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>> ES(mat);
  eigvals = ES.eigenvalues();
  return params_out;
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
parameters diagonalize(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, std::vector<norm_t<T>>& eigvals_in,
		       rokko::parameters const& params) {
  parameters params_out;
  if (get_matrix_part(params) == 'U') {
    throw std::invalid_argument("eigen3::diagonalize() : Eigen3's SelfAdjointEigenSolver does not support upper part of matrix.");
  }
  std::size_t dim = mat.rows();
  if (eigvals_in.size() < dim) eigvals_in.resize(dim);
  Eigen::Map<Eigen::Vector<norm_t<T>>> eigvals(eigvals_in.data(), eigvals_in.size());
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>> ES(mat);
  eigvals = ES.eigenvalues();
  return params_out;
}

// eigenvalues/eigenvectors
template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, int SIZE>
parameters diagonalize(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, Eigen::Vector<norm_t<T>, SIZE> & eigvals,
		       Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& eigvecs, rokko::parameters const& params) {
  parameters params_out;
  if (get_matrix_part(params) == 'U') {
    throw std::invalid_argument("eigen3::diagonalize() : Eigen3's SelfAdjointEigenSolver does not support upper part of matrix.");
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>> ES(mat);
  eigvals = ES.eigenvalues();
  eigvecs = ES.eigenvectors();
  return params_out;
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
parameters diagonalize(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, std::vector<norm_t<T>>& eigvals_in,
                Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& eigvecs, rokko::parameters const& params) {
  parameters params_out;
  if (get_matrix_part(params) == 'U') {
    throw std::invalid_argument("eigen3::diagonalize() : Eigen3's SelfAdjointEigenSolver does not support upper part of matrix.");
  }
  std::size_t dim = mat.rows();
  if (eigvals_in.size() < dim) eigvals_in.resize(dim);
  Eigen::Map<Eigen::Vector<norm_t<T>>> eigvals(eigvals_in.data(), eigvals_in.size());
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>> ES(mat);
  eigvals = ES.eigenvalues();
  eigvecs = ES.eigenvectors();
  return params_out;
}


} // namespace eigen3
} // namespace rokko

#endif // ROKKO_EIGEN3_DIAGONALIZE_HPP
