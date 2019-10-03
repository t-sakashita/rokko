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

#include <rokko/localized_matrix.hpp>
#include <rokko/dense.h>

void rokko_localized_matrix_construct(rokko_localized_matrix* matrix, int dim1, int dim2,
				      int matrix_major) {
  if (matrix_major == rokko_matrix_col_major)
    matrix->ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>(dim1, dim2);
  else
    matrix->ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(dim1, dim2);
  matrix->major = matrix_major;
}

void rokko_localized_matrix_destruct(rokko_localized_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix->ptr);
  else
    delete static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix->ptr);
  matrix->ptr = nullptr;
}

void rokko_localized_matrix_print(rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    std::cout << *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr) << std::endl;
  else
    std::cout << *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr) << std::endl;
}

void rokko_localized_matrix_generate_function(struct rokko_localized_matrix* matrix,
					      double (*func)(int i, int j)) {
  if (matrix->major == rokko_matrix_col_major)
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix->ptr), func);
  else
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix->ptr), func);
}

double* rokko_localized_matrix_get_pointer(struct rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return &(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr))(0,0);
  else
    return &(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr))(0,0);
}
