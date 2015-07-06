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
#include <rokko/eigen_exa/diagonalize_eigen_s.hpp>
#include <rokko/eigen_exa/diagonalize_eigen_sx.hpp>

namespace rokko {

namespace eigen_exa {

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
    return mapping_bc(dim_global, 1, lld, length_array, matrix_col_major_d, g);  // block_size = 1
  }

  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
		   rokko::parameters const& params, timer& timer);
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
		   VEC& eigvals,
		   parameters const& params, timer& timer);
};

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params, timer& timer) {
  if ((routine=="tri") || (routine=="eigen_s")) {
    rokko::eigen_exa::diagonalize_eigen_s(mat, eigvals, eigvecs, params, timer);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    rokko::eigen_exa::diagonalize_eigen_sx(mat, eigvals, eigvecs, params, timer);
  } else {
    std::cerr << "error: " << routine << " is not EigenExa routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
void solver::diagonalize(std::string const& routine, distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params, timer& timer) {
  if ((routine=="tri") || (routine=="eigen_s")) {
    rokko::eigen_exa::diagonalize_eigen_s(mat, eigvals, params, timer);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    rokko::eigen_exa::diagonalize_eigen_sx(mat, eigvals, params, timer);
  } else {
    std::cerr << "error: " << routine << " is not EigenExa routine" << std::endl;
    throw;
  }
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_CORE_HPP
