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

#ifndef ROKKO_UTILITY_MATRIX012_HPP
#define ROKKO_UTILITY_MATRIX012_HPP

#include <cmath>
#include <stdexcept>
#include <rokko/config.h>
#include <rokko/eigen3.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

class matrix012 {
public:

  template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
  static int get_index(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat);

  template<typename T, int ROWS, int COLS>
  static int get_index(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat, int global_i, int global_j) {
    return global_i + mat.rows() * global_j;
  }

  template<typename T, int ROWS, int COLS>
  static int get_index(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat, int global_i, int global_j) {
    return mat.cols() * global_i + global_j;
  }
  
  template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
  static void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat) {
    int n = mat.rows();
    for(int i = 0; i < mat.rows(); ++i) {
      for(int j = 0; j < mat.cols(); ++j) {
        mat(i,j) = get_index(mat, i, j);
      }
    }
  }

  template<typename T, typename MATRIX_MAJOR>
  static int get_index(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat);

  template<typename T>
  static int get_index(rokko::distributed_matrix<T, rokko::matrix_col_major>& mat, int global_i, int global_j) {
    return global_i + mat.get_m_global() * global_j;
  }

  template<typename T>
  static int get_index(rokko::distributed_matrix<T, rokko::matrix_row_major>& mat, int global_i, int global_j) {
    return mat.get_n_global() * global_i + global_j;
  }

  template<typename T, typename MATRIX_MAJOR>
  static void generate(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
    for (int local_i = 0; local_i < mat.get_m_local(); ++local_i) {
      for (int local_j = 0; local_j < mat.get_n_local(); ++local_j) {
        int global_i = mat.translate_l2g_row(local_i);
        int global_j = mat.translate_l2g_col(local_j);
        mat.set_local(local_i, local_j, get_index(mat, global_i, global_j));
      }
    }
  }
};

} // namespace rokko

#endif // ROKKO_UTILITY_MATRIX012_HPP
