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

#include <rokko/localized_vector.hpp>
#include <rokko/rokko_dense.h>

void rokko_localized_vector_construct(rokko_localized_vector* vec, int dim) {
  vec->ptr = new rokko::localized_vector<double>(dim);
}

void rokko_localized_vector_destruct(rokko_localized_vector* vec) {
  delete static_cast<rokko::localized_vector<double>*>(vec->ptr);
  vec->ptr = 0;
}

double rokko_localized_vector_get(rokko_localized_vector vec, int i) {
  return (*static_cast<rokko::localized_vector<double>*>(vec.ptr))[i];
}

/* offset by one */
double rokko_localized_vector_get_f(rokko_localized_vector vec, int i) {
  return (*static_cast<rokko::localized_vector<double>*>(vec.ptr))[i-1];
}
