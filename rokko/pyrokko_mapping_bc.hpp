/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_MAPPING_BC_HPP
#define PYROKKO_MAPPING_BC_HPP

#include <boost/variant.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/mapping_bc.hpp>
#include <rokko/pyrokko_grid.hpp>
#include <rokko/pyrokko_matrix_major_enum.hpp>
#include <rokko/utility/tuple_to_array.hpp>

namespace rokko {

class wrap_mapping_bc {
public:
  
  wrap_mapping_bc(matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _map = mapping_bc<matrix_col_major>();
    else
      _map = mapping_bc<matrix_row_major>();
  }

  wrap_mapping_bc(int global_dim, int block_size, wrap_grid const& g, matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _map = mapping_bc<matrix_col_major>(global_dim, block_size, g);
    else
      _map = mapping_bc<matrix_row_major>(global_dim, block_size, g);
  }

  wrap_mapping_bc(int global_dim, int block_size, int lld, wrap_grid const& g, matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _map = mapping_bc<matrix_col_major>(global_dim, block_size, lld, g);
    else
      _map = mapping_bc<matrix_row_major>(global_dim, block_size, lld, g);
  }

  wrap_mapping_bc(std::tuple<int,int> const& global_size, std::tuple<int,int> const& block_size, wrap_grid const& g, matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _map = mapping_bc<matrix_col_major>(to_array(global_size), to_array(block_size), g);
    else
      _map = mapping_bc<matrix_row_major>(to_array(global_size), to_array(block_size), g);
  }

  wrap_mapping_bc(mapping_bc<matrix_col_major> const& map) : is_col(true) {
    _map = mapping_bc<matrix_col_major>(map);
  }

  wrap_mapping_bc(mapping_bc<matrix_row_major> const& map) : is_col(false) {
    _map = mapping_bc<matrix_row_major>(map);
  }

  mapping_bc<matrix_col_major>& col_ver() {
    return boost::get<mapping_bc<matrix_col_major>>(_map);
  }
  
  mapping_bc<matrix_row_major>& row_ver() {
    return boost::get<mapping_bc<matrix_row_major>>(_map);
  }
  
  mapping_bc<matrix_col_major> const& col_ver() const {
    return boost::get<mapping_bc<matrix_col_major>>(_map);
  }
  
  mapping_bc<matrix_row_major> const& row_ver() const {
    return boost::get<mapping_bc<matrix_row_major>>(_map);
  }

  void set_blacs_descriptor() {
    if (is_col)
      col_ver().set_blacs_descriptor();
    else
      row_ver().set_blacs_descriptor();
  }

  bool is_col_major() const { return is_col; }

  std::string get_major_string() const {
    return is_col ? "col" : "row";
  }

  int get_mb() const {
    return is_col ? col_ver().get_mb() : row_ver().get_mb();
  }

  int get_nb() const {
    return is_col ? col_ver().get_nb() : row_ver().get_nb();
  }
  
  int get_m_global() const {
    return is_col ? col_ver().get_m_global() : row_ver().get_m_global();
  }

  int get_n_global() const {
    return is_col ? col_ver().get_n_global() : row_ver().get_n_global();
  }

  int get_m_local() const {
    return is_col ? col_ver().get_m_local() : row_ver().get_m_local();
  }

  int get_n_local() const {
    return is_col ? col_ver().get_n_local() : row_ver().get_n_local();
  }

  // local to global index
  int translate_l2g_row(int local_i) const {
    return is_col ? col_ver().translate_l2g_row(local_i) : row_ver().translate_l2g_row(local_i);
  }

  int translate_l2g_col(int local_j) const {
    return is_col ? col_ver().translate_l2g_col(local_j) : row_ver().translate_l2g_col(local_j);
  }

  // global to local index
  int translate_g2l_row(int local_i) const {
    return is_col ? col_ver().translate_g2l_row(local_i) : row_ver().translate_g2l_row(local_i);
  }

  int translate_g2l_col(int local_j) const {
    return is_col ? col_ver().translate_g2l_col(local_j) : row_ver().translate_g2l_col(local_j);
  }

  bool has_global_row_index(int global_i) const {
    return is_col ? col_ver().has_global_row_index(global_i) : row_ver().has_global_row_index(global_i);
  }

  bool has_global_col_index(int global_j) const {
    return is_col ? col_ver().has_global_col_index(global_j) : row_ver().has_global_col_index(global_j);
  }

  bool has_global_indices(std::tuple<int,int> const& global) const {
    return is_col ?
      col_ver().has_global_indices(to_array(global)) : row_ver().has_global_indices(to_array(global));
  }
  
  std::tuple<int,int> get_block_shape() const {
    return is_col ? col_ver().get_block_size() : row_ver().get_block_size();
  }
  
  std::tuple<int,int> get_global_shape() const {
    return is_col ? col_ver().get_global_size() : row_ver().get_global_size();
  }
  
  std::tuple<int,int> get_local_shape() const {
    return is_col ? col_ver().get_local_size() : row_ver().get_local_size();
  }

  std::tuple<int,int> translate_l2g(std::tuple<int,int> const& local) const {
    return is_col ? col_ver().translate_l2g(to_array(local)) : row_ver().translate_l2g(to_array(local));
  }

  std::tuple<int,int> translate_g2l(std::tuple<int,int> const& global) const {
    return is_col ? col_ver().translate_g2l(to_array(global)) : row_ver().translate_g2l(to_array(global));
  }

  wrap_grid get_grid() const {
    const grid& g = is_col ? col_ver().get_grid() : row_ver().get_grid();
    return wrap_grid(g);
  }

private:
  bool is_col;
  boost::variant<mapping_bc<matrix_row_major>, mapping_bc<matrix_col_major>> _map;
};

} // end namespace rokko

#endif // PYROKKO_MAPPING_BC_HPP
