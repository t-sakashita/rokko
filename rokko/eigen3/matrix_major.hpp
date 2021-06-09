/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/matrix_major.hpp>
#include <Eigen/Dense>

namespace rokko {

namespace detail {

template<typename T>
struct major_t {};

template<typename T, int ROWS>
struct major_t<Eigen::Matrix<T, ROWS, 1, Eigen::ColMajor>> {
  using major_type = rokko::matrix_col_major;
};

template<typename T, int ROWS, int COLS>
struct major_t<Eigen::Matrix<T, ROWS, COLS, Eigen::ColMajor>> {
  using major_type = rokko::matrix_col_major;
};

template<typename T, int ROWS, int COLS>
struct major_t<Eigen::Matrix<T, ROWS, COLS, Eigen::RowMajor>> {
  using major_type = rokko::matrix_row_major;
};

} // end namespace detail

template<typename MATRIX_MAJOR>
using matrix_major = typename detail::major_t<MATRIX_MAJOR>::major_type;

namespace detail {

template<typename MATRIX_MAJOR>
struct eigen3_matrix_major;

template<>
struct eigen3_matrix_major<rokko::matrix_row_major> {
  static constexpr int value = Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major> {
  static constexpr int value = Eigen::ColMajor;
};

} // end namespace detail

template<typename MATRIX_MAJOR>
constexpr int eigen3_major = detail::eigen3_matrix_major<MATRIX_MAJOR>::value;

} // namespace rokko
