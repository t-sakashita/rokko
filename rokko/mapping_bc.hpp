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

template <typename MATRIX_MAJOR>
class mapping_bc : public mapping_global2local, public mapping_local2array<MATRIX_MAJOR> {
public:
  typedef MATRIX_MAJOR matrix_major;
  explicit mapping_bc() {}
  explicit mapping_bc(int global_dim, int block_size)
    : mapping_global2local(global_dim, block_size, grid()), mapping_local2array<MATRIX_MAJOR>() {}
  explicit mapping_bc(int global_dim, int block_size, MATRIX_MAJOR matrix_major, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array<MATRIX_MAJOR>() {}
  // This consutuctor is used by default_mapping func of solver interface for elpa, elemental, scalapack
  
  explicit mapping_bc(int global_dim, int block_size, int padded_row_size, int padded_col_size, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array<MATRIX_MAJOR>(padded_row_size, padded_col_size) {}
  // This consutuctor is used by default_mapping func of solver interface for eigen_exa
  
  explicit mapping_bc(int global_dim, int block_size, int lld, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array<MATRIX_MAJOR>(lld) {}
  // This consutuctor can be used for every solvers
  
  explicit mapping_bc(int global_dim, int block_size, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in), mapping_local2array<MATRIX_MAJOR>() {}
  // defalut major is col-major
  
  //explicit mapping_bc(int m_global_in, int n_global_in, int mb_in, int nb_in, grid const& g_in)
  //  : mapping_global2local(m_global_in, n_global_in, mb_in, nb_in, g_in), mapping_local2array<MATRIX_MAJOR>() {}

  template<typename SOLVER>
  explicit mapping_bc(int global_dim, grid const& g_in, SOLVER const& solver_in) {
    //std::cout << "constructor: m_local=" << solver_in.default_mapping(g_in, global_dim).get_m_local() << " n_local=" << solver_in.default_mapping(g_in, global_dim).get_n_local() << std::endl;
    *this = solver_in.default_mapping(global_dim, g_in);
  }

};

} // namespace rokko

#endif // ROKKO_MAPPING_BC_HPP
