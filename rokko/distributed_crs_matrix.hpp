/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/anasazi/mapping_1d.hpp>

namespace rokko {

namespace detail {

class distributed_crs_matrix_base {
public:
  virtual void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) = 0;
  virtual void insert(int row, int col_size, int* cols, double* const values) = 0;
  virtual void complete() = 0;
  virtual int get_dim() = 0;
  virtual int num_local_rows() = 0;
  virtual int start_row() = 0;
  virtual int end_row() = 0;
  virtual void print() = 0;
};

} // end namespace detail

class distributed_crs_matrix {
public:
  template<typename SOLVER>
  explicit distributed_crs_matrix(int row_dim, int col_dim, SOLVER& solver_in) {
    mat = solver_in.create_distributed_crs_matrix(row_dim, col_dim);
  }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    mat->insert(row, cols, values);
  }
  void insert(int row, int col_size, int* cols, double* const values) {
    mat->insert(row, col_size, cols, values);
  }
  void complete() {
    mat->complete();
  }
  int get_dim() {
    return mat->get_dim();
  }
  int num_local_rows() {
    return mat->num_local_rows();
  }
  int start_row() {
    return mat->start_row();
  }
  int end_row() {
    return mat->end_row();
  }
  void print() {
    mat->print();
  }
  detail::distributed_crs_matrix_base* get_matrix() {
    return mat;
  }

private:
  detail::distributed_crs_matrix_base* mat;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

