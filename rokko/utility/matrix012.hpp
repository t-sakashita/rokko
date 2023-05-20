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

  template<int MATRIX_MAJOR>
  static int get_index(int lld, int global_i, int global_j);

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat) {
    const auto ld = mat.rows();

    for(int j = 0; j < mat.cols(); ++j) {
      for(int i = 0; i < mat.rows(); ++i) {
        mat(i,j) = get_index<Eigen::ColMajor>(ld, i, j);
      }
    }
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat) {
    const auto ld = mat.cols();

    for(int i = 0; i < mat.rows(); ++i) {
      for(int j = 0; j < mat.cols(); ++j) {
        mat(i,j) = get_index<Eigen::RowMajor>(ld, i, j);
      }
    }
  }

  static int get_index(rokko::mapping_bc<rokko::matrix_col_major> const& map, int global_i, int global_j) {
    return global_i + map.get_m_global() * global_j;
  }

  static int get_index(rokko::mapping_bc<rokko::matrix_row_major> const& map, int global_i, int global_j) {
    return map.get_n_global() * global_i + global_j;
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_col_major>& mat) {
    const auto& map = mat.get_mapping();

    for (int local_j = 0; local_j < map.get_n_local(); ++local_j) {
      int global_j = map.translate_l2g_col(local_j);
      for (int local_i = 0; local_i < map.get_m_local(); ++local_i) {
        int global_i = map.translate_l2g_row(local_i);
        mat.set_local(local_i, local_j, get_index(map, global_i, global_j));
      }
    }
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_row_major>& mat) {
    const auto& map = mat.get_mapping();

    for (int local_i = 0; local_i < map.get_m_local(); ++local_i) {
      int global_i = map.translate_l2g_row(local_i);
      for (int local_j = 0; local_j < map.get_n_local(); ++local_j) {
        int global_j = map.translate_l2g_col(local_j);
        mat.set_local(local_i, local_j, get_index(map, global_i, global_j));
      }
    }
  }
};

template<>
int matrix012::get_index<Eigen::ColMajor>(int ld, int global_i, int global_j) {
  return global_i + ld * global_j;
}

template<>
int matrix012::get_index<Eigen::RowMajor>(int ld, int global_i, int global_j) {
  return ld * global_i + global_j;
}

} // namespace rokko
