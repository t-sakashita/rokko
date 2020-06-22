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

#ifndef ROKKO_ELPA_DIAG_HPP
#define ROKKO_ELPA_DIAG_HPP

#include <rokko/elpa/elpa.h>

namespace rokko {
namespace elpa {

int diag(elpa_t handle, double* a, double* ev, double* q) {
  int error;
  elpa_eigenvectors_d(handle, a, ev, q, &error);
  return error;
}

int diag(elpa_t handle, float* a, float* ev, float* q) {
  int error;
  elpa_eigenvectors_f(handle, a, ev, q, &error);
  return error;
}

int diag(elpa_t handle, double* a, double* ev) {
  int error;
  elpa_eigenvalues_d(handle, a, ev, &error);
  return error;
}

int diag(elpa_t handle, float* a, float* ev) {
  int error;
  elpa_eigenvalues_f(handle, a, ev, &error);
  return error;
}

template<typename MATRIX, typename VECTOR>
int diag(elpa_t handle, MATRIX& a, VECTOR& ev, MATRIX& q) {
  return diag(handle, a.get_array_pointer(), storage(ev), q.get_array_pointer());
}

template<typename MATRIX, typename VECTOR>
int diag(elpa_t handle, MATRIX& a, VECTOR& ev) {
  return diag(handle, a.get_array_pointer(), storage(ev));
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAG_HPP
