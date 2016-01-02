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

#ifndef ROKKO_ELPA_CORE_HPP
#define ROKKO_ELPA_CORE_HPP

#include <rokko/parameters.hpp>
#include <rokko/elpa/diagonalize_elpa1.hpp>
#include <rokko/elpa/diagonalize_elpa2.hpp>

namespace rokko {
namespace elpa {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return true; }
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const {
    // Determine mb, nb, lld, larray
    int mb = dim / g.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = dim / g.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int b = std::min(mb, nb);
    return mapping_bc<matrix_col_major>(dim, b, g);
  }
  template <typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat, VEC& eigvals,
			 distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params) {
    std::string routine = "";
    if(params.defined("routine")) {
      routine = params.get_string("routine");
    }
    if (routine=="elpa1") {
      return rokko::elpa::diagonalize_elpa1(mat, eigvals, eigvecs, params);
    } else if (routine=="elpa2") {
      return rokko::elpa::diagonalize_elpa2(mat, eigvals, eigvecs, params);
    } else if (routine=="") {  // default
      return rokko::elpa::diagonalize_elpa1(mat, eigvals, eigvecs, params);
    } else {
      std::cerr << "error: " << routine << " is not ELPA's routine" << std::endl;
      throw;
    }
  }
  template <typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat, VEC& eigvals,
			 parameters const& params) {
    std::cerr << "not yet implemented" << std::endl;
    throw;
    //return rokko::elpa::diagonalize(mat, eigvals, params);
  }
};

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_CORE_HPP
