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

#ifndef PYROKKO_DISTRIBUTED_MATRIX_HPP
#define PYROKKO_DISTRIBUTED_MATRIX_HPP

#include <boost/variant.hpp>
#include <pybind11/pybind11.h>

#include <rokko/pyrokko_mapping_bc.hpp>
#include <rokko/distributed_matrix.hpp>


namespace rokko {

class wrap_distributed_matrix {
public:

  wrap_distributed_matrix(wrap_mapping_bc const& map_in) : is_col(map_in.is_col_major()), map(map_in) {
    if (is_col)
      _ptr = new distributed_matrix<double,matrix_col_major>(map_in.col_ver());
    else
      _ptr = new distributed_matrix<double,matrix_row_major>(map_in.row_ver());
  }

  distributed_matrix<double,matrix_col_major>& col_ver() {
    return *(boost::get<distributed_matrix<double,matrix_col_major>*>(_ptr));
  }
  
  distributed_matrix<double,matrix_row_major>& row_ver() {
    return *(boost::get<distributed_matrix<double,matrix_row_major>*>(_ptr));
  }

  distributed_matrix<double,matrix_col_major> const& col_ver() const {
    return *(boost::get<distributed_matrix<double,matrix_col_major>*>(_ptr));
  }
  
  distributed_matrix<double,matrix_row_major> const& row_ver() const {
    return *(boost::get<distributed_matrix<double,matrix_row_major>*>(_ptr));
  }

  // map member function
  int get_mb() const {
    return is_col ? col_ver().get_mb() : row_ver().get_mb();
  }

  int get_nb() const {
    return is_col ? col_ver().get_nb() : row_ver().get_nb();
  }
  
  py::tuple get_block_shape() const {
    return py::make_tuple(get_mb(), get_nb());
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
  int translate_l2g_row(const int& local_i) const {
    return is_col ? col_ver().translate_l2g_row(local_i) : row_ver().translate_l2g_row(local_i);
  }

  int translate_l2g_col(const int& local_j) const {
    return is_col ? col_ver().translate_l2g_col(local_j) : row_ver().translate_l2g_col(local_j);
  }

  // global to local index
  int translate_g2l_row(const int& local_i) const {
    return is_col ? col_ver().translate_g2l_row(local_i) : row_ver().translate_g2l_row(local_i);
  }

  int translate_g2l_col(const int& local_j) const {
    return is_col ? col_ver().translate_g2l_col(local_j) : row_ver().translate_g2l_col(local_j);
  }

  bool is_gindex_myrow(const int& global_i) const {
    return is_col ? col_ver().is_gindex_myrow(global_i) : row_ver().is_gindex_myrow(global_i);
  }

  bool is_gindex_mycol(const int& global_j) const {
    return is_col ? col_ver().is_gindex_mycol(global_j) : row_ver().is_gindex_mycol(global_j);
  }

  bool is_gindex(const int& global_i, const int& global_j) const {
    return is_col ?
      col_ver().is_gindex(global_i, global_j) : row_ver().is_gindex(global_i, global_j);
  }
  
  py::tuple get_global_shape() const {
    return py::make_tuple(get_m_global(), get_n_global());
  }
  
  py::tuple get_local_shape() const {
    return py::make_tuple(get_m_local(), get_n_local());
  }

  std::tuple<int,int> translate_l2g(std::tuple<int,int> const& local) const {
    return std::make_tuple(translate_l2g_row(std::get<0>(local)),
                           translate_l2g_col(std::get<1>(local)));
  }

  std::tuple<int,int> translate_g2l(std::tuple<int,int> const& global) const {
    return std::make_tuple(translate_g2l_row(std::get<0>(global)),
                           translate_g2l_col(std::get<1>(global)));
  }
  
  void set_local(int local_i, int local_j, double value) {
    if (is_col)
      col_ver().set_local(local_i, local_j, value);
    else
      row_ver().set_local(local_i, local_j, value);
  }
  
  void set_global(int global_i, int global_j, double value) {
    if (is_col)
      col_ver().set_global(global_i, global_j, value);
    else
      row_ver().set_global(global_i, global_j, value);
  }

  int get_length_array() const {
    return is_col ? col_ver().get_length_array() : row_ver().get_length_array();
  }
  int get_lld() const {
    return is_col ? col_ver().get_lld() : row_ver().get_lld();
  };
  
  void set_zeros() {
    if (is_col)
      col_ver().set_zeros();
    else
      row_ver().set_zeros();
  }

  void generate(std::function<double(int, int)> const& func) {
    if (is_col)
      col_ver().generate(func);
    else
      row_ver().generate(func);
  }
  
  void print() const {
    if (is_col)
      col_ver().print();
    else
      row_ver().print();
  }

  bool is_major_col() const {
    return is_col;
  }
  
  std::string get_major_string() const {
    return is_col ? "col" : "row";
  }

  wrap_mapping_bc& get_map() {
    return map;
  }

private:
  bool is_col;
  boost::variant<distributed_matrix<double,matrix_row_major>*, distributed_matrix<double,matrix_col_major>*> _ptr;
  wrap_mapping_bc map;
};

} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MATRIX_HPP
