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
#include <rokko/parameters.hpp>
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

  mapping_bc<matrix_col_major> optimized_mapping(int global_dim, grid const& g) const {
    int nx, ny;
    int nprow = g.get_nprow();
    int npcol = g.get_npcol();
    int n = global_dim;
    ROKKO_eigen_exa_get_matdims(nprow, npcol, n, &nx, &ny);
    //std::cout << "nx=" << nx << std::endl;
    int lld = nx;
    return mapping_bc<matrix_col_major>(global_dim, 1, lld, g);  // block_size = 1
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

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="tri") || (routine=="eigen_s")) {
    return rokko::eigen_exa::diagonalize_eigen_s(mat, eigvals, eigvecs, params);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    return rokko::eigen_exa::diagonalize_eigen_sx(mat, eigvals, eigvecs, params);
  } else {
    std::cerr << "error: " << routine << " is not EigenExa routine" << std::endl;
    throw;
  }
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  std::string routine = "";
  if(params.defined("routine")) {
    routine = params.get_string("routine");
  }
  if ((routine=="tri") || (routine=="eigen_s")) {
    return rokko::eigen_exa::diagonalize_eigen_s(mat, eigvals, params);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    return rokko::eigen_exa::diagonalize_eigen_sx(mat, eigvals, params);
  } else {
    std::cerr << "error: " << routine << " is not EigenExa routine" << std::endl;
    throw;
  }
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_CORE_HPP
