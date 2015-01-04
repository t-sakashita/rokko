/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN_EXA_CORE_HPP
#define ROKKO_EIGEN_EXA_CORE_HPP

#include <rokko/eigen_exa/eigen_exa.h>
#include <rokko/eigen_exa/diagonalize.hpp>

namespace rokko {
namespace eigen_exa {

struct eigen_s {};
struct eigen_sx {};

template<typename ROUTINE>
class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return true; }
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  template<typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    int n = mat.get_m_global();

    // calculate sizes of my proc's local part of distributed matrix
    int NPROW = mat.get_grid().get_nprow();
    int NPCOL = mat.get_grid().get_npcol();

    int n1 = ((n-1)/NPROW+1);
    int nm;
    int i1 = 6, i2 = 16*4, i3 = 16*4*2;
    ROKKO_cstab_get_optdim(n1, i1, i2, i3, &nm);

    int NB  = 64;
    int nmz = ((n-1)/NPROW+1);
    nmz = ((nmz-1)/NB+1)*NB+1;
    int nmw = ((n-1)/NPCOL+1);
    nmw = ((nmw-1)/NB+1)*NB+1;

    int larray = std::max(nmz, nm)*nmw;
    
    mat.set_lld(nm);
    mat.set_length_array(larray);
    mat.set_block_size(1, 1);
    int m_local = mat.calculate_row_size();
    int n_local = mat.calculate_col_size();
    mat.set_local_size(m_local, n_local);
  }
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, VEC& eigvals,
		   distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer);
};

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::eigen_exa::eigen_s>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::eigen_exa::diagonalize_s(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::eigen_exa::eigen_sx>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
    rokko::eigen_exa::diagonalize_sx(mat, eigvals, eigvecs, timer);
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_CORE_HPP
