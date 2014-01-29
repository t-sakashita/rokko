/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/localized_vector.hpp>
#include <rokko/rokko.h>

void rokko_localized_vector_construct(rokko_localized_vector* vec, int dim) {
  vec->ptr = new rokko::localized_vector(dim);
}

void rokko_localized_vector_destruct(rokko_localized_vector* vec) {
  delete static_cast<rokko::localized_vector*>(vec->ptr);
}

double rokko_localized_vector_get(rokko_localized_vector vec, int i) {
  return (*static_cast<rokko::localized_vector*>(vec.ptr))[i];
}

/* offset by one */
double rokko_localized_vector_get_f(rokko_localized_vector vec, int i) {
  return (*static_cast<rokko::localized_vector*>(vec.ptr))[i-1];
}
