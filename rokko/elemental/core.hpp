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
#include <rokko/parameters.hpp>
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

  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const {
    return mapping_bc<matrix_col_major>(dim, 1, matrix_col_major_d, g);  // block_size = 1
  }

  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params) {
    std::string routine = "";
    if(params.defined("routine")) {
      routine = params.get_string("routine");
    }
    if ((routine == "") || (routine == "pmrrr")) {
      return rokko::elemental::diagonalize(mat, eigvals, eigvecs, params);
    } else {
      throw std::invalid_argument("elemental::diagonalize() : " + routine + " is invalid routine name");
    }
  }
  
  template<typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
			 VEC& eigvals,
			 parameters const& params) {
    std::string routine = "";
    if(params.defined("routine")) {
      routine = params.get_string("routine");
    }
    if ((routine == "") || (routine == "pmrrr")) {
      return rokko::elemental::diagonalize(mat, eigvals, params);
    } else {
      throw std::invalid_argument("elemental::diagonalize() : " + routine + " is invalid routine name");
    }
  }
};

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_CORE_HPP
