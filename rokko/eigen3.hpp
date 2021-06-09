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

#pragma once

#include <rokko/eigen3/matrix_major.hpp>
#include <Eigen/Dense>

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
int ld(Eigen::Matrix<T, ROWS, COLS, MATRIX_MAJOR> const& mat) {
  return mat.innerSize();
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
