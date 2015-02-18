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

#ifndef ROKKO_ELEMENTAL_CORE_HPP
#define ROKKO_ELEMENTAL_CORE_HPP

#include <El.hpp>
#include <rokko/elemental/diagonalize.hpp>
#include <boost/type_traits/is_same.hpp>

namespace rokko {
namespace elemental {

class solver {
public:
  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) {
    return boost::is_same<GRID_MAJOR, grid_col_major_t>::value;
  }
  void initialize(int& argc, char**& argv) { El::Initialize(argc, argv); }
  void finalize() { El::Finalize(); }
  // void optimized_grid_size() {}

  mapping_bc optimized_mapping(grid const& g, int dim) const {
    return mapping_bc(g, dim, 1);  // block_size = 1
  }

  template<typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    mat.set_default_lld();
    mat.set_default_length_array();
  }

  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, VEC& eigvals,
                   distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
    rokko::elemental::diagonalize(mat, eigvals, eigvecs, timer);
  }
};

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_CORE_HPP
