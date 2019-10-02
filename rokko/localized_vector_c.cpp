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

void rokko_localized_vector_construct(rokko_localized_vector* vec, int dim) {
  vec->ptr = new Eigen::VectorXd(dim);
}

void rokko_localized_vector_destruct(rokko_localized_vector* vec) {
  delete static_cast<Eigen::VectorXd*>(vec->ptr);
  vec->ptr = nullptr;
}

double rokko_localized_vector_get(rokko_localized_vector vec, int i) {
  return (*static_cast<Eigen::VectorXd*>(vec.ptr))[i];
}

/* offset by one */
double rokko_localized_vector_get_f(rokko_localized_vector vec, int i) {
  return (*static_cast<Eigen::VectorXd*>(vec.ptr))[i-1];
}

void rokko_localized_vector_print(rokko_localized_vector vec) {
  std::cout << *static_cast<Eigen::VectorXd*>(vec.ptr) << std::endl;
}
