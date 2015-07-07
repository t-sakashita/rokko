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
  explicit mapping_bc(int global_dim, int block_size, MATRIX_MAJOR matrix_major)
    : mapping_global2local(global_dim, block_size, grid()), mapping_local2array(matrix_major) {}
  template<typename MATRIX_MAJOR>
  explicit mapping_bc(int global_dim, int block_size, MATRIX_MAJOR matrix_major, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array(matrix_major) {}
  // This consutuctor is used by optimized_mapping func of solver interface for elpa, elemental, scalapack

  template<typename MATRIX_MAJOR>
  explicit mapping_bc(int global_dim, int block_size, int lld, MATRIX_MAJOR matrix_major, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array(lld, matrix_major) {}
  // This consutuctor is used by optimized_mapping func of solver interface for eigen_exa
  
  template<typename MATRIX_MAJOR>
  explicit mapping_bc(int m_global_in, int n_global_in, int mb_in, int nb_in, MATRIX_MAJOR matrix_major, grid const& g_in)
    : mapping_global2local(m_global_in, n_global_in, mb_in, nb_in, g_in), mapping_local2array(matrix_major) {}

  explicit mapping_bc(int global_dim, int block_size, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array(matrix_col_major_d) {}
  // defalut major is col-major

  template<typename SOLVER>
  explicit mapping_bc(int global_dim, grid const& g_in, SOLVER const& solver_in) {
    //std::cout << "constructor: m_local=" << solver_in.optimized_mapping(g_in, global_dim).get_m_local() << " n_local=" << solver_in.optimized_mapping(g_in, global_dim).get_n_local() << std::endl;
    *this = solver_in.optimized_mapping(global_dim, g_in);
  }

private:
};

} // namespace rokko

#endif // ROKKO_MAPPING_BC_HPP
