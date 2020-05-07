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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_QR_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_QR_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/cscalapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/scalapack/diagonalize_psyevx.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

// qr (pdsyevx) eigenvalues & eigenvectors
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_qr(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			       parameters params) {
  if (params.defined("abstol")) {
    T abstol = params.get<T>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  parameters params_out = diagonalize_psyevx(mat, eigvals, eigvecs, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished pdsyevx (qr)" << std::endl;
  }

  return params_out;
}

// qr (pdsyevx) only eigenvalues
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_qr(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters params) {
  if (params.defined("abstol")) {
    T abstol = params.get<T>("abstol");
    if (abstol > 0) {
      params.set("abstol", - abstol);
    }
  }
  parameters params_out = diagonalize_psyevx(mat, eigvals, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished pdsyevx (qr)" << std::endl;
  }

  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_QR_HPP
