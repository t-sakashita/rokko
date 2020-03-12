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

#ifndef ROKKO_EIGENEXA_SOLVER_HPP
#define ROKKO_EIGENEXA_SOLVER_HPP

#include <rokko/matrix_major.hpp>
#include <rokko/parameters.hpp>
#include <rokko/eigenexa/diagonalize_eigen_s.hpp>
#include <rokko/eigenexa/diagonalize_eigen_sx.hpp>

namespace rokko {

namespace eigenexa {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& /* grid_major */) { return true; }
  void initialize(int& /* argc */, char**& /* argv */) {}
  void finalize() {}

  mapping_bc<matrix_col_major> default_mapping(int global_dim, grid const& g) const {
    int nx, ny;
    std::tie(nx, ny) = rokko::eigenexa::get_matdims(g, global_dim);
    return mapping_bc<matrix_col_major>(global_dim, 1, {nx, ny}, g);  // block_size = 1
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
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="tri") || (routine=="eigen_s")) {
    return rokko::eigenexa::diagonalize_eigen_s(mat, eigvals, eigvecs, params);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    return rokko::eigenexa::diagonalize_eigen_sx(mat, eigvals, eigvecs, params);
  } else {
    throw std::invalid_argument("eigenexa::diagonalize() : " + routine + " is invalid routine name");
  }
}

template<typename MATRIX_MAJOR, typename VEC>
parameters solver::diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if ((routine=="tri") || (routine=="eigen_s")) {
    return rokko::eigenexa::diagonalize_eigen_s(mat, eigvals, params);
  } else if ((routine=="") || (routine=="penta") || (routine=="eigen_sx")) {
    return rokko::eigenexa::diagonalize_eigen_sx(mat, eigvals, params);
  } else {
    throw std::invalid_argument("eigenexa::diagonalize() : " + routine + " is invalid routine name");
  }
}

} // namespace eigenexa
} // namespace rokko

#endif // ROKKO_EIGENEXA_SOLVER_HPP
