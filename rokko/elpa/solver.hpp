/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/parameters.hpp>
#include <rokko/default_block_size.hpp>
#include <rokko/elpa/elpa.h>
#include <rokko/elpa/diagonalize.hpp>

namespace rokko {
namespace elpa {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& /* grid_major */) { return true; }
  void initialize(int& /* argc */, char**& /* argv */) {
    if (elpa_init(20200417) != ELPA_OK) {
      throw std::invalid_argument("ERROR: elpa::initialize()");
    }
  }
  void finalize() {
    // int error;
    // elpa_uninit(&error);
  }
  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const {
    // Determine mb, nb, lld, larray
    const auto mb = get_default_block_size(dim, g.get_nprow());
    const auto nb = get_default_block_size(dim, g.get_npcol());
    // Note: it should be that mb = nb in pdsyev.
    const auto b = std::min(mb, nb);
    return mapping_bc<matrix_col_major>(dim, b, g);
  }
  template <typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat, VEC& eigvals,
			 distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			 parameters const& params) {
    return rokko::elpa::diagonalize(mat, eigvals, eigvecs, params);
  }

  template <typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat, VEC& eigvals,
			 parameters const& params) {
    return rokko::elpa::diagonalize(mat, eigvals, params);
  }
};

} // namespace elpa
} // namespace rokko
