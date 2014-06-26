/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
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
#include <iostream>

namespace rokko {
namespace scalapack {

struct pdsyev {};
struct pdsyevd {};
struct pdsyevx {};

template<typename ROUTINE>
class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return true; }
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  void optimized_grid_size() {}
  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    // Determine mb, nb, lld, larray
    int mb = mat.get_m_global() / mat.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = mat.get_n_global() / mat.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int tmp = std::min(mb, nb);
    mat.set_block_size(tmp, tmp);

    // Determine m_local, n_local from m_global, n_global, mb, nb
    mat.set_default_local_size();
    mat.set_default_lld();
    mat.set_default_length_array();
  }
  template<typename MATRIX_MAJOR, typename TIMER>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   distributed_matrix<MATRIX_MAJOR>& eigvecs, TIMER& timer);
};

template<>
template<typename MATRIX_MAJOR, typename TIMER>
void solver<rokko::scalapack::pdsyev>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat,
  localized_vector& eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, TIMER& timer) {
  rokko::scalapack::diagonalize(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename TIMER>
void solver<rokko::scalapack::pdsyevd>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat,
  localized_vector& eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, TIMER& timer) {
  rokko::scalapack::diagonalize_d(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename TIMER>
void solver<rokko::scalapack::pdsyevx>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat,
  localized_vector& eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, TIMER& timer) {
  rokko::scalapack::diagonalize_x(mat, eigvals, eigvecs, timer);
}

} // namespace sclapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_CORE_HPP
