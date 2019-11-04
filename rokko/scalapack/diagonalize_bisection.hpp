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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_BISECTION_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_BISECTION_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/cscalapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/scalapack/diagonalize_pdsyevx.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

// bisection (pdsyevx) eigenvalues / eigenvectors
template<typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_bisection(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<double>("abstol") < 0) {
      std::stringstream msg;
      msg << "scalapack::diagonalize_bisection() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use pdsyevx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }
  parameters params_out = diagonalize_pdsyevx(mat, eigvals, eigvecs, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished pdsyevx (bisection)" << std::endl;
  }

  return params_out;
}

// bisection (pdsyevx) only eigenvalues
template<typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_bisection(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  if (params.defined("abstol")) {
    if (params.get<double>("abstol") < 0) {
      std::stringstream msg;
      msg << "scalapack::diagonalize_bisection() : " << std::endl
          << "abstol is negative value, which means QR method." << std::endl
          << "To use pdsyevx as bisection solver, set abstol a positive value" << std::endl;
      throw std::invalid_argument(msg.str());
    }
  }
  parameters params_out = diagonalize_pdsyevx(mat, eigvals, params);

  if (params.get_bool("verbose")) {
    std::cout << "finished pdsyevx (bisection)" << std::endl;
  }

  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_BISECTION_HPP
