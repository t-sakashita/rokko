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

#include <rokko/elpa/elpa.hpp>
#include <rokko/elpa/diagonalize.hpp>

namespace rokko {
namespace elpa {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) { return true; }
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  //void optimized_grid_size() {}
  mapping_bc optimized_mapping(grid const& g, int dim)  const {
    // Determine mb, nb, lld, larray
    int mb = dim / g.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = dim / g.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int b = std::min(mb, nb);
    return mapping_bc(g, dim, b);
  }
  template<typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<double, MATRIX_MAJOR>& mat) {
    mat.set_default_lld();
    mat.set_default_length_array();
  }
  template <typename MATRIX_MAJOR, typename VEC>
  void diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat, VEC& eigvals,
                   distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
    rokko::elpa::diagonalize(mat, eigvals, eigvecs, timer);
  }
};

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_CORE_HPP
