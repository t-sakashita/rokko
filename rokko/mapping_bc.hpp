/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>,
*                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_BC_HPP
#define ROKKO_MAPPING_BC_HPP

#include <rokko/grid.hpp>
#include <rokko/matrix_major.hpp>
#include <rokko/mapping_global2local.hpp>
#include <rokko/mapping_local2array.hpp>

namespace rokko {

class mapping_bc : public mapping_global2local, public mapping_local2array {
public:
  explicit mapping_bc() {}
  template<typename MATRIX_MAJOR>
  explicit mapping_bc(grid const& g_in, int global_dim, int block_size, MATRIX_MAJOR matrix_major)
    : mapping_global2local(g_in, global_dim, global_dim, block_size), mapping_local2array(matrix_major) {}

  template<typename MATRIX_MAJOR>
  explicit mapping_bc(grid const& g_in, int global_dim, int block_size, int lld, int length_array, MATRIX_MAJOR matrix_major)
    : mapping_local2array(lld, length_array, matrix_major), mapping_global2local(g_in, global_dim, global_dim, block_size) {}
  
  template<typename MATRIX_MAJOR>
  explicit mapping_bc(grid const& g_in, int m_global_in, int n_global_in, int block_size, MATRIX_MAJOR matrix_major)
    : mapping_global2local(g_in, m_global_in, n_global_in, block_size), mapping_local2array(matrix_major) {}

  explicit mapping_bc(grid const& g_in, int global_dim, int block_size)
    : mapping_global2local(g_in, global_dim, global_dim, block_size), mapping_local2array(matrix_col_major_d) {}  // default of matrix major is matrix_col_major

  template<typename SOLVER>
  explicit mapping_bc(grid const& g_in, int global_dim, SOLVER const& solver_in) {
    //std::cout << "constructor: m_local=" << solver_in.optimized_mapping(g_in, global_dim).get_m_local() << " n_local=" << solver_in.optimized_mapping(g_in, global_dim).get_n_local() << std::endl;

    *this = solver_in.optimized_mapping(g_in, global_dim);
  }

private:
};

} // namespace rokko

#endif // ROKKO_MAPPING_BC_HPP
