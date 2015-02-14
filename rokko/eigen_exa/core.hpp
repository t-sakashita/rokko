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
  int get_lld(int NPROW, int NPCOL,  {

    }
  template<typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    int nx, ny;
    int dim = mat.get_mapping().get_dim();
    int NPROW = mat.get_mapping().get_nprow();
    int NPCOL = mat.get_mapping().get_npcol();  
    n1 = ((n-1)/NPROW+1);  
    int i1 = 6, i2 = 16*4, i3 = 16*4*2; // input for cstab_get_optdim
    int nm; // output for cstab_get_optdim
    ROKKO_cstab_get_optdim(n1, i1, i2, i3, &nm);

    NB  = eigen_NB;

    nmz = ((n-1)/NPROW+1);
    nmz = ((nmz-1)/NB+1)*NB+1;
                 nn  = nmz
    nmz = (n-1)/NB+1
    nmz = ((nmz-1)/NPROW+1)*NB
    nmz = MAX(nn, nmz)

    nmw = ((n-1)/NPCOL+1)
    nmw = ((nmw-1)/NB+1)*NB+1
                 nn  = nmw
    nmw = (n-1)/NB+1
    nmw = ((nmw-1)/NPCOL+1)*NB
    nmw = MAX(nn, nmw)

    int larray = MAX(nmz, nm)*nmw

    int nx = nm;
    int ny = (larray-1)/nm+1;
    //ROKKO_eigen_get_matdims( dim, &nx, &ny );

    std::cout << "nx=" << nx << std::endl;
    std::cout << "larray=" << larray << std::endl;
    mat.set_lld(nx);
    mat.set_length_array(larray);
  }
  mapping_bc optimized_mapping(grid const& g, int dim) const {
    return mapping_bc(g, dim, 1);  // block_size = 1
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
