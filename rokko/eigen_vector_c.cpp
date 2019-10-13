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

#include <rokko/eigen3.hpp>
#include <rokko/dense.h>

void rokko_eigen_vector_construct(rokko_eigen_vector* vec, int dim) {
  vec->ptr = new Eigen::VectorXd(dim);
}

void rokko_eigen_vector_construct_array_sizes(struct rokko_eigen_vector* vector,
                                              int dim, double* ptr) {
  vector->ptr = new Eigen::Map<Eigen::VectorXd>(ptr, dim);
}

void rokko_eigen_vector_destruct(rokko_eigen_vector* vec) {
  delete static_cast<Eigen::VectorXd*>(vec->ptr);
  vec->ptr = nullptr;
}

int rokko_eigen_vector_get_dim(rokko_eigen_vector vec) {
  return static_cast<Eigen::VectorXd*>(vec.ptr)->size();
}

double* rokko_eigen_vector_get_array_pointer(rokko_eigen_vector vec) {
  return static_cast<Eigen::VectorXd*>(vec.ptr)->data();
}

double rokko_eigen_vector_get(rokko_eigen_vector vec, int i) {
  return (*static_cast<Eigen::VectorXd*>(vec.ptr))[i];
}

/* offset by one */
double rokko_eigen_vector_get_f(rokko_eigen_vector vec, int i) {
  return (*static_cast<Eigen::VectorXd*>(vec.ptr))[i-1];
}

void rokko_eigen_vector_print(rokko_eigen_vector vec) {
  std::cout << *static_cast<Eigen::VectorXd*>(vec.ptr) << std::endl;
}
