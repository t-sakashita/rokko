/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/collective.hpp>
#include <rokko/dense.h>

void rokko_gather(struct rokko_distributed_matrix matrix, double* array, int root) {
  if (matrix.major == rokko_matrix_col_major){
    rokko::distributed_matrix<double, rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr);
    rokko::gather(*ptr_, array, root);
  } else {
    rokko::distributed_matrix<double, rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr);
    rokko::gather(*ptr_, array, root);
  }
}

void rokko_scatter(double* global_array, struct rokko_distributed_matrix matrix, int root) {
  if (matrix.major == rokko_matrix_col_major){
    rokko::distributed_matrix<double, rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr);
    rokko::scatter(global_array, *ptr_, root);
  } else {
    rokko::distributed_matrix<double, rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr);
    rokko::scatter(global_array, *ptr_, root);
  }
}

void rokko_gather_localized_matrix(struct rokko_distributed_matrix matrix, struct rokko_localized_matrix lmatrix, int root) {
  if (matrix.major == rokko_matrix_col_major){
    rokko::distributed_matrix<double, rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>* loc_ptr_ = static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(lmatrix.ptr);
    rokko::gather(*ptr_, *loc_ptr_, root);
  } else {
    rokko::distributed_matrix<double, rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>* loc_ptr_ = static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(lmatrix.ptr);
    rokko::gather(*ptr_, *loc_ptr_, root);
  }
}

void rokko_scatter_localized_matrix(struct rokko_localized_matrix lmatrix, struct rokko_distributed_matrix matrix, int root) {
  if (matrix.major == rokko_matrix_col_major){
    rokko::distributed_matrix<double, rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>* loc_ptr_ = static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(lmatrix.ptr);
    rokko::scatter(*loc_ptr_, *ptr_, root);
  } else {
    rokko::distributed_matrix<double, rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>* loc_ptr_ = static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(lmatrix.ptr);
    rokko::scatter(*loc_ptr_, *ptr_, root);
  }
}

void rokko_all_gather(struct rokko_distributed_matrix matrix, double* array) {
  if (matrix.major == rokko_matrix_col_major){
    rokko::distributed_matrix<double, rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr);
    for(int root = 0; root < ptr_->get_nprocs(); ++root) {
      rokko::gather(*ptr_, array, root);
    }
  } else {
    rokko::distributed_matrix<double, rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr);
    for(int root = 0; root < ptr_->get_nprocs(); ++root) {
      rokko::gather(*ptr_, array, root);
    }
  }
}

