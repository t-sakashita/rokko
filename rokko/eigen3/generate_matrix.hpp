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

#ifndef ROKKO_EIGEN3_GENERATE_MATRIX_HPP
#define ROKKO_EIGEN3_GENERATE_MATRIX_HPP

#include <iostream>
#include <rokko/eigen3.hpp>

namespace rokko {

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

#endif // ROKKO_EIGEN3_GENERATE_MATRIX_HPP
