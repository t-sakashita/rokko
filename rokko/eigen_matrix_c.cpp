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

#include <rokko/eigen_matrix.h>
#include <rokko/eigen3.hpp>
#include <rokko/eigen3/generate_matrix.hpp>

void rokko_eigen_matrix_construct(rokko_eigen_matrix* matrix, int dim1, int dim2,
				      int matrix_major) {
  if (matrix_major == rokko_matrix_col_major)
    matrix->ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>(dim1, dim2);
  else
    matrix->ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(dim1, dim2);
  matrix->major = matrix_major;
}

void rokko_eigen_matrix_construct_array_sizes(struct rokko_eigen_matrix* matrix,
                                              int dim1, int dim2, double* ptr, int matrix_major) {
  if (matrix_major == rokko_matrix_col_major) {
    matrix->ptr = new Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(ptr, dim1, dim2);
  } else {
    matrix->ptr = new Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(ptr, dim1, dim2);
  }
  matrix->major = matrix_major;
}

void rokko_eigen_matrix_destruct(rokko_eigen_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix->ptr);
  else
    delete static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix->ptr);
  matrix->ptr = nullptr;
}

void rokko_eigen_matrix_print(rokko_eigen_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    std::cout << *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr) << std::endl;
  else
    std::cout << *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr) << std::endl;
}

void rokko_eigen_matrix_generate_function(struct rokko_eigen_matrix matrix,
					      double (*func)(int i, int j)) {
  if (matrix.major == rokko_matrix_col_major)
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr), func);
  else
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr), func);
}

void rokko_eigen_matrix_generate_function_p(struct rokko_eigen_matrix matrix,
					      double (*func)(const int* i, const int* j)) {

  const auto g = [&func](int i, int j) { return func(&i, &j); };

  if (matrix.major == rokko_matrix_col_major)
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr), g);
  else
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr), g);
}

void rokko_eigen_matrix_generate_function1(struct rokko_eigen_matrix matrix,
					      double (*func)(int i, int j)) {
  const auto g = [&func](int i, int j) { return func(i+1, j+1); };

  if (matrix.major == rokko_matrix_col_major)
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr), g);
  else
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr), g);
}

void rokko_eigen_matrix_generate_function1_p(struct rokko_eigen_matrix matrix,
					      double (*func)(const int* i, const int* j)) {
  const auto g = [&func](int i, int j) {
    const int i1 = i+1, j1 = j+1;
    return func(&i1, &j1);
  };

  if (matrix.major == rokko_matrix_col_major)
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr), g);
  else
    rokko::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr), g);
}

double rokko_eigen_matrix_get(struct rokko_eigen_matrix matrix, int i, int j) {
  if (matrix.major == rokko_matrix_col_major)
    return (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr))(i, j);
  else
    return (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr))(i, j);
}

double rokko_eigen_matrix_get1(struct rokko_eigen_matrix matrix, int i, int j) {
  if (matrix.major == rokko_matrix_col_major)
    return (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr))(i-1, j-1);
  else
    return (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr))(i-1, j-1);
}

void rokko_eigen_matrix_set(struct rokko_eigen_matrix matrix, int i, int j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr))(i, j) = value;
  else
    (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr))(i, j) = value;
}

void rokko_eigen_matrix_set1(struct rokko_eigen_matrix matrix, int i, int j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr))(i-1, j-1) = value;
  else
    (*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr))(i-1, j-1) = value;
}

int rokko_eigen_matrix_get_m(struct rokko_eigen_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr)->rows();
  else
    return static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr)->rows();
}

int rokko_eigen_matrix_get_n(struct rokko_eigen_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr)->cols();
  else
    return static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr)->cols();
}

bool rokko_eigen_matrix_is_row_major(struct rokko_eigen_matrix matrix) {
  return matrix.major == rokko_matrix_row_major;
}

bool rokko_eigen_matrix_is_col_major(struct rokko_eigen_matrix matrix) {
  return matrix.major == rokko_matrix_col_major;
}

double* rokko_eigen_matrix_get_array_pointer(struct rokko_eigen_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return storage(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr));
  else
    return storage(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr));
}
