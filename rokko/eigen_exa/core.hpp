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

#ifndef ROKKO_EIGEN_EXA_CORE_HPP
#define ROKKO_EIGEN_EXA_CORE_HPP

#include <rokko/matrix_major.hpp>
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

  mapping_bc optimized_mapping(grid const& g, int dim_global) const {
    int nx, ny;
    int nprow = g.get_nprow();
    int npcol = g.get_npcol();
    int n = dim_global;
    ROKKO_eigen_exa_get_matdims(nprow, npcol, n, &nx, &ny);

    //std::cout << "nx=" << nx << std::endl;
    int lld = nx;
    int length_array = nx * ny;
    return mapping_bc(g, dim_global, 1, lld, length_array, matrix_col_major_d);  // block_size = 1
  }

  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat, VEC& eigvals,
		   distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer);
};

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::eigen_exa::eigen_s>::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::eigen_exa::diagonalize_s(mat, eigvals, eigvecs, timer);
}

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::eigen_exa::eigen_sx>::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
  VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
    rokko::eigen_exa::diagonalize_sx(mat, eigvals, eigvecs, timer);
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_CORE_HPP
