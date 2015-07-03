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

#include <rokko/scalapack/diagonalize.hpp>
#include <rokko/scalapack/diagonalize_pdsyevd.hpp>
#include <rokko/scalapack/diagonalize_pdsyevx.hpp>

namespace rokko {
namespace scalapack {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return true; }
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  mapping_bc optimized_mapping(grid const& g, int dim)  const {
    // Determine mb, nb, lld, larray
    int mb = dim / g.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = dim / g.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int b = std::min(mb, nb);
    return mapping_bc(dim, b, g);
  }

  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
		   rokko::parameters const& params, timer& timer);
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
		   rokko::parameters const& params, timer& timer);
};

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params, timer& timer) {
  /*if ((routine=="dsyev") || (routine=="qr")) {
    rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params, timer);
  } else if ((routine=="dsyevr") || (routine=="mr3")) {
    rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params, timer);
  } else if ((routine=="dsyevd") || (routine=="dc")) {
    rokko::lapack::diagonalize_dsyevd(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="dsyevx") {
    rokko::lapack::diagonalize_dsyevx(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="bisection") {
    rokko::lapack::diagonalize_bisection(mat, eigvals, eigvecs, params, timer);
  } else if (routine=="") {
    if (params.defined("upper_limit") || params.defined("lower_limit")) {
      rokko::lapack::diagonalize_dsyevr(mat, eigvals, eigvecs, params, timer);
    }
    else {
      rokko::lapack::diagonalize_dsyev(mat, eigvals, eigvecs, params, timer);
    }
  } else {
    std::cerr << "error: " << routine << " is not scalapack routine" << std::endl;
    throw;
  }*/
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 rokko::parameters const& params, timer& timer) {
}

/*
template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::scalapack::pdsyev>::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::scalapack::diagonalize(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::scalapack::pdsyevd>::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::scalapack::diagonalize_d(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::scalapack::pdsyevx>::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::scalapack::diagonalize_x(mat, eigvals, eigvecs, timer);
}
*/

} // namespace sclapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_CORE_HPP
