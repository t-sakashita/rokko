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
#include <rokko/eigen3/generate_vector.hpp>
#include <rokko/dense.h>

void rokko_eigen_vector_construct(rokko_eigen_vector* vec, int dim) {
  vec->ptr = new Eigen::VectorXd(dim);
}

void rokko_eigen_vector_construct_array_size(struct rokko_eigen_vector* vector,
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

void rokko_eigen_vector_set(rokko_eigen_vector vec, int i, double val) {
  (*static_cast<Eigen::VectorXd*>(vec.ptr))[i] = val;
}

/* offset by one */
void rokko_eigen_vector_set_f(rokko_eigen_vector vec, int i, double val) {
  (*static_cast<Eigen::VectorXd*>(vec.ptr))[i-1] = val;
}

void rokko_eigen_vector_print(rokko_eigen_vector vec) {
  std::cout << *static_cast<Eigen::VectorXd*>(vec.ptr) << std::endl;
}

void rokko_eigen_vector_generate_function(rokko_eigen_vector vec,
					      double (*func)(int i)) {
  rokko::generate(*static_cast<Eigen::VectorXd*>(vec.ptr), func);
}

void rokko_eigen_vector_generate_function_p(rokko_eigen_vector vec,
                                            double (*func)(const int* i)) {
  rokko::generate(*static_cast<Eigen::VectorXd*>(vec.ptr),
                  [&func](int i) { return func(&i); } );
}

void rokko_eigen_vector_generate_function_f(rokko_eigen_vector vec,
                                            double (*func)(int i)) {
  rokko::generate(*static_cast<Eigen::VectorXd*>(vec.ptr),
                  [&func](int i) { return func(i+1); });
}

void rokko_eigen_vector_generate_function_f_p(rokko_eigen_vector vec,
                                              double (*func)(const int* i)) {
  rokko::generate(*static_cast<Eigen::VectorXd*>(vec.ptr),
                  [&func](int i) {
                    int i1 = i+1;
                    return func(&i1);
                  } );
}
