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

#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/parameters.hpp>
#include <rokko/scalapack/diagonalize_pdsyev.hpp>
#include <rokko/scalapack/diagonalize_pdsyevx.hpp>
#include <rokko/scalapack/diagonalize_pdsyevd.hpp>
#ifdef ROKKO_HAVE_PDSYEVR
#include <rokko/scalapack/diagonalize_pdsyevr.hpp>
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

  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params);
  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params);
};

// eigenvalues / eigenvectors
template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="pdsyev") || (routine=="qr")) {
    return rokko::scalapack::diagonalize_pdsyev(mat, eigvals, eigvecs, params);
  } else if ((routine=="pdsyevr") || (routine=="mr3")) {
#ifdef ROKKO_HAVE_PDSYEVR
    return rokko::scalapack::diagonalize_pdsyevr(mat, eigvals, eigvecs, params);
#else
    throw std::invalid_argument("scalapack::diagonalize() : the routine pdsyevr does not exist in your machine.");
#endif
  } else if ((routine=="pdsyevd") || (routine=="dc")) {
    return rokko::scalapack::diagonalize_pdsyevd(mat, eigvals, eigvecs, params);
  } else if (routine=="pdsyevx") {
    return rokko::scalapack::diagonalize_pdsyevx(mat, eigvals, eigvecs, params);
    //} else if (routine=="bisection") {
    //rokko::scalapack::diagonalize_bisection(mat, eigvals, eigvecs, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
#ifdef ROKKO_HAVE_PDSYEVR
      return rokko::scalapack::diagonalize_pdsyevr(mat, eigvals, eigvecs, params);
#else
      throw std::invalid_argument("scalapack::diagonalize() : the default routine for a range of eigenvalues, pdsyevr does not exist in your machine.");
#endif
    } else {
      return rokko::scalapack::diagonalize_pdsyev(mat, eigvals, eigvecs, params);
    }
  } else {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " is invalid routine name");
  }
}

// only eigenvalues
template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="pdsyev") || (routine=="qr")) {
    return rokko::scalapack::diagonalize_pdsyev(mat, eigvals, params);
  } else if ((routine=="pdsyevr") || (routine=="mr3")) {
#ifdef ROKKO_HAVE_PDSYEVR
    return rokko::scalapack::diagonalize_pdsyevr(mat, eigvals, params);
#else
    throw std::invalid_argument("scalapack::diagonalize() : the routine pdsyevr does not exist in your machine.");
#endif
  } else if ((routine=="pdsyevd") || (routine=="dc")) {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " does not support computing only eigenvalues");
  } else if (routine=="pdsyevx") {
    return rokko::scalapack::diagonalize_pdsyevx(mat, eigvals, params);
    //} else if (routine=="bisection") {
    //rokko::scalapack::diagonalize_bisection(mat, eigvals, params);
  } else if (routine=="") {
    if (lapack::is_interval(params)) {
#ifdef ROKKO_HAVE_PDSYEVR
      return rokko::scalapack::diagonalize_pdsyevr(mat, eigvals, params);
#else
      throw std::invalid_argument("scalapack::diagonalize() : the default routine for a range of eigenvalues, pdsyevr does not exist in your machine.");
#endif
    } else {
      return rokko::scalapack::diagonalize_pdsyev(mat, eigvals, params);
    }
  } else {
    throw std::invalid_argument("scalapack::diagonalize() : " + routine + " is invalid routine name");
  }
}

} // namespace sclapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_CORE_HPP
