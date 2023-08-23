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

#pragma once

#include <rokko/elpa/elpa.h>


namespace rokko {
namespace elpa {

int deallocate(const elpa_t handle) {
  int error;
  elpa_deallocate(handle, &error);
  return error;
}

int set(const elpa_t handle, const char *name, int value) {
  int error;
  elpa_set_integer(handle, name, value, &error);
  return error;
}

int set(const elpa_t handle, const char *name, double value) {
  int error;
  elpa_set_double(handle, name, value, &error);
  return error;
}


int diag(const elpa_t handle, double* a, double* ev, double* q) {
  int error;
  elpa_eigenvectors(handle, a, ev, q, &error);
  return error;
}

int diag(const elpa_t handle, float* a, float* ev, float* q) {
  int error;
  elpa_eigenvectors(handle, a, ev, q, &error);
  return error;
}

int diag(const elpa_t handle, std::complex<double>* a, double* ev, std::complex<double>* q) {
  int error;
  elpa_eigenvectors(handle, a, ev, q, &error);
  return error;
}

int diag(const elpa_t handle, std::complex<float>* a, float* ev, std::complex<float>* q) {
  int error;
  elpa_eigenvectors(handle, a, ev, q, &error);
  return error;
}

int diag(const elpa_t handle, double* a, double* ev) {
  int error;
  elpa_eigenvalues(handle, a, ev, &error);
  return error;
}

int diag(const elpa_t handle, float* a, float* ev) {
  int error;
  elpa_eigenvalues(handle, a, ev, &error);
  return error;
}

int diag(const elpa_t handle, std::complex<double>* a, double* ev) {
  int error;
  elpa_eigenvalues(handle, a, ev, &error);
  return error;
}

int diag(const elpa_t handle, std::complex<float>* a, float* ev) {
  int error;
  elpa_eigenvalues(handle, a, ev, &error);
  return error;
}

template<typename MATRIX, typename VECTOR>
int diag(const elpa_t handle, MATRIX& a, VECTOR& ev, MATRIX& q) {
  std::cout << "AAAAAAA\n";
  return diag(handle, a.get_array_pointer(), storage(ev), q.get_array_pointer());
}

template<typename MATRIX, typename VECTOR>
int diag(const elpa_t handle, MATRIX& a, VECTOR& ev) {
  return diag(handle, a.get_array_pointer(), storage(ev));
}

} // namespace elpa
} // namespace rokko
