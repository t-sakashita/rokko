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
#include <stdexcept>
#include <rokko/eigen3.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

// tridiagonal_toeplitz_matrix can be defined for n >= 2.
class tridiagonal_toeplitz_matrix {
public:
  template<typename T>
  static void multiply(int dim, const T* x, T* y, T a, T b) {
    y[0] = a * x[0] + b * x[1];
    y[dim-1] = a * x[dim-1] + b * x[dim - 2];
    for (int k = 1; k < (dim-1); ++k) { // from 1 to end-1
      y[k] = b * x[k-1] + a * x[k] + b * x[k+1];
    }
  }

  template<typename T>
  static void multiply(int dim, const std::vector<T>& v, std::vector<T>& w, T a, T b) {
    multiply(dim, v.data(), w.data());
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat, T a, T b) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("tridiagonal_toeplitz_matrix::generate() : non-square matrix");
    mat.setZero();
    const int n = mat.cols();
    mat(0, 0) = a;  mat(1, 0) = b;
    for(int i = 1; i < n-1; ++i) {
      mat(i-1, i) = b;
      mat(i,   i) = a;
      mat(i+1, i) = b;
    }
    mat(n-2, n-1) = b;  mat(n-1, n-1) = a;
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat, T a, T b) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("tridiagonal_toeplitz_matrix::generate() : non-square matrix");
    mat.setZero();
    const int n = mat.rows();
    mat(0, 0) = a;  mat(0, 1) = b;
    for(int i = 1; i < n-1; ++i) {
      mat(i, i-1) = b;
      mat(i, i) = a;
      mat(i, i+1) = b;
    }
    mat(n-1, n-2) = b;  mat(n-1, n-1) = a;
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_col_major>& mat, T a, T b) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("tridiagonal_toeplitz_matrix::generate() : non-square matrix");
    mat.set_zeros();
    const int n = mat.get_n_global();
    mat.set_global(0, 0, a);  mat.set_global(1, 0, b);
    for(int i = 1; i < n-1; ++i) {
      mat.set_global(i-1, i, b);
      mat.set_global(i,   i, a);
      mat.set_global(i+1, i, b);
    }
    mat.set_global(n-2, n-1, b);  mat.set_global(n-1, n-1, a);
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_row_major>& mat, T a, T b) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("tridiagonal_toeplitz_matrix::generate() : non-square matrix");
    mat.set_zeros();
    const int n = mat.get_m_global();
    mat.set_global(0, 0, a);  mat.set_global(0, 1, b);
    for(int i = 1; i < n-1; ++i) {
      mat.set_global(i, i-1, b);
      mat.set_global(i, i, a);
      mat.set_global(i, i+1, b);
    }
    mat.set_global(n-1, n-2, b);  mat.set_global(n-1, n-1, a);
  }

  // calculate k-th smallest eigenvalue of dim-dimensional Tridiagonal_Toeplitz matrix (k=0...dim-1)
  template<typename T>
  static T eigenvalue(int dim, int k, T a, T b) {
    return a - 2 * b * std::cos(M_PI * (k + 1) / (dim + 1));
  }
};

} // namespace rokko
