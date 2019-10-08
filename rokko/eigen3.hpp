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

#ifndef ROKKO_EIGEN3_HPP
#define ROKKO_EIGEN3_HPP

#include <rokko/matrix_major.hpp>
#include <rokko/traits/major_t.hpp>
#include <Eigen/Dense>

namespace rokko {

template<typename T, int ROWS>
struct major_t<Eigen::Matrix<T, ROWS, 1, Eigen::ColMajor>> {
  typedef rokko::matrix_col_major major_type;
};

template<typename T, int ROWS, int COLS>
struct major_t<Eigen::Matrix<T, ROWS, COLS, Eigen::ColMajor>> {
  typedef rokko::matrix_col_major major_type;
};

template<typename T, int ROWS, int COLS>
struct major_t<Eigen::Matrix<T, ROWS, COLS, Eigen::RowMajor>> {
  typedef rokko::matrix_row_major major_type;
};

} // namespace rokko

namespace Eigen {

template<typename T, int ROWS, int MATRIX_MAJOR>
int size(Eigen::Matrix<T, ROWS, 1, MATRIX_MAJOR> const& vec) {
  return vec.size();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
bool is_col_major(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const&) {
  return MATRIX_MAJOR == Eigen::ColMajor;
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
bool is_row_major(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const&) {
  return MATRIX_MAJOR == Eigen::RowMajor;
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
int rows(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return mat.rows();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
int cols(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return mat.cols();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
int lda(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return is_col_major(mat) ? rows(mat) : cols(mat);
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
const T* storage(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return mat.data();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
T* storage(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR>& mat) {
  return mat.data();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
const T* storage(Eigen::Ref<Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR>> const& mat) {
  return mat.data();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
T* storage(Eigen::Ref<Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR>>& mat) {
  return mat.data();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
const std::complex<T>* storage(Eigen::Matrix<std::complex<T>, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return mat.data();
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
std::complex<T>* storage(Eigen::Matrix<std::complex<T>, ROWS, COLS, MATRIX_MAJOR>& mat) {
  return mat.data();
}

template<typename T, int ROWS = Eigen::Dynamic>
using Vector = Eigen::Matrix<T, ROWS, 1>;

template<typename T>
using RefVec = Eigen::Ref<Vector<T>>;

} // namespace Eigen

#endif // ROKKO_EIGEN3_HPP
