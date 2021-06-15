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

#pragma once

#include <iostream>
#include <rokko/eigen3.hpp>

namespace rokko {

template<typename T, int ROWS, int COLS, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat, FUNC func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

template<typename T, int ROWS, int COLS, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat, FUNC func) {
  for(int j = 0; j < mat.cols(); ++j) {
    for(int i = 0; i < mat.rows(); ++i) {
      mat(i, j) = func(i, j);
    }
  }
}

template<typename T, int ROWS, int COLS>
void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat, std::function<T(int, int)> const& func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

template<typename T, int ROWS, int COLS>
void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat, std::function<T(int, int)> const& func) {
  for(int j = 0; j < mat.cols(); ++j) {
    for(int i = 0; i < mat.rows(); ++i) {
      mat(i, j) = func(i, j);
    }
  }
}

} // namespace rokko
