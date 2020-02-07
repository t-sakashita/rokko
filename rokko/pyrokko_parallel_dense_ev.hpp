/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_PARALLEL_DENSE_EV_HPP
#define PYROKKO_PARALLEL_DENSE_EV_HPP

#include <rokko/pyrokko_mapping_bc.hpp>
#include <rokko/pyrokko_distributed_matrix.hpp>
#include <rokko/parallel_dense_ev.hpp>

namespace rokko {

class wrap_parallel_dense_ev : public parallel_dense_ev {
public:
  wrap_parallel_dense_ev(std::string const& solver_name) : parallel_dense_ev(solver_name) {}
  
  wrap_parallel_dense_ev() = default;

  void initialize() {
    int num = 1;
    char** ptr = NULL;
    parallel_dense_ev::initialize(num, ptr);
  }

  wrap_mapping_bc default_mapping(int dim, wrap_grid const& g) const {
    return wrap_mapping_bc(parallel_dense_ev::default_mapping(dim, g.get_grid()));
  }

  template<typename VEC>
  wrap_parameters diagonalize(wrap_distributed_matrix& mat,	VEC& eigvals, wrap_distributed_matrix& eigvecs,
			 wrap_parameters const& params) {
    assert(mat.is_major_col() == eigvecs.is_major_col());
    if (mat.is_major_col())
      return parallel_dense_ev::diagonalize(mat.col_ver(), eigvals, eigvecs.col_ver(), parameters(params));
    else
      return parallel_dense_ev::diagonalize(mat.row_ver(), eigvals, eigvecs.row_ver(), parameters(params));
  }

  template<typename VEC>
  wrap_parameters diagonalize(wrap_distributed_matrix& mat,	VEC& eigvals,
			 wrap_parameters const& params) {
    if (mat.is_major_col())
      return parallel_dense_ev::diagonalize(mat.col_ver(), eigvals, parameters(params));
    else
      return parallel_dense_ev::diagonalize(mat.row_ver(), eigvals, parameters(params));
  }

};

} // end namespace rokko

#endif // PYROKKO_PARALLEL_DENSE_EV_HPP
