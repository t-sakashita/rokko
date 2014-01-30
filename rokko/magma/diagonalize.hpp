/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAGMA_DIAGONALIZE_H
#define ROKKO_MAGMA_DIAGONALIZE_H

#include "rokko/localized_matrix.hpp"
#include "rokko/localized_vector.hpp"

// Use CUBLAS for now
#include <rokko/config.hpp>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>  // include before magma.h
#include <magma.h>
#include <magma_lapack.h>

namespace rokko {
namespace magma {

template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, double* eigvals, localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int info;

  int dim = mat.rows();
  magma_int_t lwork, liwork, aux_iwork[1];
  double aux_work[1];

  magma_dsyevd(MagmaVec, 'V', dim, NULL, dim, NULL, aux_work, -1, aux_iwork, -1, &info);
  lwork  = (magma_int_t) aux_work[0];
  double h_work[lwork];
  liwork = aux_iwork[0];

  magma_int_t *iwork;
  // eigenvalue decomposition
  timer_in.start(1);
  magma_dsyevd(MagmaVec, 'V', dim, &mat(0,0), dim, eigvals, h_work, lwork, iwork, liwork, &info);
  timer_in.stop(1);
  eigvecs = mat;

  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }

  return info;
}

template<class MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace magma
} // namespace rokko

#endif // ROKKO_MAGMA_DIAGONALIZE_H
