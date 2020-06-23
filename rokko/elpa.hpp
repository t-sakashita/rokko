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
#include <complex>
#undef complex

namespace rokko {
namespace elpa {

int elpa_set(elpa_t handle, const char *name, int value) {
  int error;
  elpa_set_integer(handle, name, value, &error);
  return error;
}

int elpa_set(elpa_t handle, const char *name, double value) {
  int error;
  elpa_set_double(handle, name, value, &error);
  return error;
}


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

int diag(elpa_t handle, std::complex<double>* a, double* ev, std::complex<double>* q) {
  int error;
  elpa_eigenvectors_dc(handle, complex_cast(a), ev, complex_cast(q), &error);
  return error;
}

int diag(elpa_t handle, std::complex<float>* a, float* ev, std::complex<float>* q) {
  int error;
  elpa_eigenvectors_fc(handle, complex_cast(a), ev, complex_cast(q), &error);
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

int diag(elpa_t handle, std::complex<double>* a, double* ev) {
  int error;
  elpa_eigenvalues_dc(handle, complex_cast(a), ev, &error);
  return error;
}

int diag(elpa_t handle, std::complex<float>* a, float* ev) {
  int error;
  elpa_eigenvalues_fc(handle, complex_cast(a), ev, &error);
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
