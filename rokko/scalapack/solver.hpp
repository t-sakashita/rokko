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

#pragma once

#include <rokko/config.h>
#include <rokko/parameters.hpp>
#include <rokko/scalapack/diagonalize_psyev.hpp>
#include <rokko/scalapack/diagonalize_psyevx.hpp>
#include <rokko/scalapack/diagonalize_psyevd.hpp>
#include <rokko/scalapack/diagonalize_qr.hpp>
#include <rokko/scalapack/diagonalize_bisection.hpp>
#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
#include <rokko/scalapack/diagonalize_psyevr.hpp>
#endif

namespace rokko {
namespace scalapack {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& /* grid_major */) { return true; }
  void initialize(int& /* argc */, char**& /* argv */) {}
  void finalize() {}
  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g)  const {
    // Determine mb, nb, lld, larray
    int mb = dim / g.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = dim / g.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int b = std::min(mb, nb);
    return mapping_bc<matrix_col_major>(dim, b, g);
  }

  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params);
};

// eigenvalues / eigenvectors
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="syev") || (routine=="qr")) {
    return rokko::scalapack::diagonalize_psyev(mat, eigvals, eigvecs, params);
  } else if ((routine=="syevr") || (routine=="mr3")) {
#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
    return rokko::scalapack::diagonalize_psyevr(mat, eigvals, eigvecs, params);
#else
    throw std::invalid_argument("scalapack::diagonalize() : the routine pdsyevr does not exist in your machine.");
#endif
  } else if ((routine=="syevd") || (routine=="dc")) {
    return rokko::scalapack::diagonalize_psyevd(mat, eigvals, eigvecs, params);
  } else if (routine=="syevx") {
    return rokko::scalapack::diagonalize_psyevx(mat, eigvals, eigvecs, params);
  } else if (routine=="bisection") {
    return rokko::scalapack::diagonalize_bisection(mat, eigvals, eigvecs, params);
  } else if (routine=="qr") {
    return rokko::scalapack::diagonalize_qr(mat, eigvals, eigvecs, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
      return rokko::scalapack::diagonalize_psyevr(mat, eigvals, eigvecs, params);
#else
      throw std::invalid_argument("scalapack::diagonalize() : the default routine for a range of eigenvalues, pdsyevr does not exist in your machine.");
#endif
    } else {
      return rokko::scalapack::diagonalize_psyev(mat, eigvals, eigvecs, params);
    }
  } else {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " is invalid routine name");
  }
}

// only eigenvalues
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="syev") || (routine=="qr")) {
    return rokko::scalapack::diagonalize_psyev(mat, eigvals, params);
  } else if ((routine=="syevr") || (routine=="mr3")) {
#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
    return rokko::scalapack::diagonalize_psyevr(mat, eigvals, params);
#else
    throw std::invalid_argument("scalapack::diagonalize() : the routine pdsyevr does not exist in your machine.");
#endif
  } else if ((routine=="syevd") || (routine=="dc")) {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " does not support computing only eigenvalues");
  } else if (routine=="syevx") {
    return rokko::scalapack::diagonalize_psyevx(mat, eigvals, params);
  } else if (routine=="bisection") {
    return rokko::scalapack::diagonalize_bisection(mat, eigvals, params);
  } else if (routine=="qr") {
    return rokko::scalapack::diagonalize_qr(mat, eigvals, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
#ifdef ROKKO_HAVE_SCALAPACK_PDSYEVR
      return rokko::scalapack::diagonalize_psyevr(mat, eigvals, params);
#else
      throw std::invalid_argument("scalapack::diagonalize() : the default routine for a range of eigenvalues, pdsyevr does not exist in your machine.");
#endif
    } else {
      return rokko::scalapack::diagonalize_psyev(mat, eigvals, params);
    }
  } else {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " is invalid routine name");
  }
}

} // namespace sclapack
} // namespace rokko
