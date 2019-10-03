/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LOCALIZED_MATRIX_HPP
#define ROKKO_LOCALIZED_MATRIX_HPP

#include <iostream>
#include <rokko/eigen3.hpp>

namespace rokko {

namespace detail {
    
template<typename MATRIX_MAJOR>
struct eigen3_matrix_major;

template<>
struct eigen3_matrix_major<rokko::matrix_row_major> {
  static const int value = Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major> {
  static const int value = Eigen::ColMajor;
};

} // end namespace detail

template<typename MATRIX_MAJOR>
constexpr int eigen3_major = detail::eigen3_matrix_major<MATRIX_MAJOR>::value;


template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, FUNC func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, std::function<T(int, int)> const& func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
