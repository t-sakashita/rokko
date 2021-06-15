/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <rokko/eigen3.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

// laplacian_matrix can be defined for n >= 2.
class laplacian_matrix {
public:
  template<typename T>
  static void multiply(int dim, const T* x, T* y) {
    y[0] = x[0] - x[1];
    y[dim-1] = 2 * x[dim-1] - x[dim - 2];
    for (int k = 1; k < (dim-1); ++k) { // from 1 to end-1
      y[k] = - x[k-1] + 2 * x[k] - x[k+1];
    }
  }

  template<typename T>
  static void multiply(int dim, const std::vector<T>& v, std::vector<T>& w) {
    multiply(dim, v.data(), w.data());
  }

  // Giving diagonal and sub-diagonal elements as two vectors
  template<typename VEC>
  static void generate(VEC& d, VEC& e) {
    if (d.size() > (e.size()+1))
      throw std::invalid_argument("laplacian_matrix::generate() : The vector for sub-diagonal elements is too short.");

    const int n = d.size();
    d[0] = 1;
    for(int i = 1; i < n; ++i)
      d[i] = 2;
    for(int i = 0; i < (n-1); ++i)
      e[i] = -1;
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("laplacian_matrix::generate() : non-square matrix");
    mat.setZero();
    const int n = mat.cols();
    mat(0, 0) = 1;  mat(1, 0) = -1;
    for(int i = 1; i < n-1; ++i) {
      mat(i-1, i) = -1;
      mat(i,   i) = 2;
      mat(i+1, i) = -1;
    }
    mat(n-2, n-1) = -1;  mat(n-1, n-1) = 2;
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("laplacian_matrix::generate() : non-square matrix");
    mat.setZero();
    const int n = mat.rows();
    mat(0, 0) = 1;  mat(0, 1) = -1;
    for(int i = 1; i < n-1; ++i) {
      mat(i, i-1) = -1;
      mat(i, i) = 2;
      mat(i, i+1) = -1;
    }
    mat(n-1, n-2) = -1;  mat(n-1, n-1) = 2;
  }

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_col_major>& mat) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("laplacian_matrix::generate() : non-square matrix");
    mat.set_zeros();
    const int n = mat.get_n_global();
    mat.set_global(0, 0, 1);  mat.set_global(1, 0, -1);
    for(int i = 1; i < n-1; ++i) {
      mat.set_global(i-1, i, -1);
      mat.set_global(i,   i, 2);
      mat.set_global(i+1, i, -1);
    }
    mat.set_global(n-2, n-1, -1);  mat.set_global(n-1, n-1, 2);
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_row_major>& mat) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("laplacian_matrix::generate() : non-square matrix");
    mat.set_zeros();
    const int n = mat.get_m_global();
    mat.set_global(0, 0, 1);  mat.set_global(0, 1, -1);
    for(int i = 1; i < n-1; ++i) {
      mat.set_global(i, i-1, -1);
      mat.set_global(i, i, 2);
      mat.set_global(i, i+1, -1);
    }
    mat.set_global(n-1, n-2, -1);  mat.set_global(n-1, n-1, 2);
  }
#endif

  // calculate k-th smallest eigenvalue of dim-dimensional Laplacian matrix (k=0...dim-1)
  static double eigenvalue(int dim, int k) {
    return 2 * (1 - std::cos(M_PI * (2 * k + 1) / (2 * dim + 1)));
  }
};

} // namespace rokko
