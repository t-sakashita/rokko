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
#include <rokko/utility/tuple_to_array.hpp>


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
  
  std::tuple<int,int> get_block_shape() const {
    return is_col ? col_ver().get_block_size() : row_ver().get_block_size();
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

  bool is_gindex_myrow(int global_i) const {
    return is_col ? col_ver().is_gindex_myrow(global_i) : row_ver().is_gindex_myrow(global_i);
  }

  bool is_gindex_mycol(int global_j) const {
    return is_col ? col_ver().is_gindex_mycol(global_j) : row_ver().is_gindex_mycol(global_j);
  }

  bool is_gindex(int global_i, int global_j) const {
    return is_col ?
      col_ver().is_gindex(global_i, global_j) : row_ver().is_gindex(global_i, global_j);
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

  template <typename MATRIX_MAJOR>
  distributed_matrix<double,MATRIX_MAJOR>* get_ptr() {
    return boost::get<distributed_matrix<double,MATRIX_MAJOR>*>(_ptr);
  }

  template <typename MATRIX_MAJOR>
  double* get_array_pointer() {
    return get_ptr<MATRIX_MAJOR>()->get_array_pointer();
  }

  py::array_t<double> get_ndarray() {
    if (is_col)
      return py::array_t<double>({get_m_local(), get_n_local()}, {sizeof(double), sizeof(double)*get_lld()}, get_array_pointer<matrix_col_major>(), py::cast(*this));
    else
      return py::array_t<double>({get_m_local(), get_n_local()}, {sizeof(double)*get_lld(), sizeof(double)}, get_array_pointer<matrix_row_major>(), py::cast(*this));
  }

  template <typename MATRIX_MAJOR>
  auto get_eigen_map() {
    return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, rokko::eigen3_major<MATRIX_MAJOR>>,0,Eigen::OuterStride<>>(get_array_pointer<MATRIX_MAJOR>(), get_m_local(), get_n_local(), Eigen::OuterStride<Eigen::Dynamic>(get_lld()));
  }

  void set_ndarray(py::array_t<double> const& mat) {
    py::array_t<double> array;
    if (is_col) {
      array = py::array_t<double>({get_m_local(), get_n_local()}, {sizeof(double), sizeof(double)*get_lld()}, get_array_pointer<matrix_col_major>(), py::cast(*this));
    } else {
      array = py::array_t<double>({get_m_local(), get_n_local()}, {sizeof(double)*get_lld(), sizeof(double)}, get_array_pointer<matrix_row_major>(), py::cast(*this));
    }

    auto r = array.template mutable_unchecked<2>();
    for (auto i = 0; i < r.shape(0); ++i) {
      for (auto j = 0; j < r.shape(1); ++j) {
        r(i, j) = *mat.data(i, j);
      }
    }
  }

  void set_col_major_matrix(Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> mat) {
    if (!is_col)
      throw std::invalid_argument("Cannot set col-major ndarray to row-major distributed_matrix");

    get_eigen_map<matrix_col_major>() = static_cast<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(mat);
  }

  void set_row_major_matrix(Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mat) {
    if (is_col)
      throw std::invalid_argument("Cannot set row-major ndarray to col-major distributed_matrix");

    get_eigen_map<matrix_row_major>() = static_cast<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mat);
  }

private:
  bool is_col;
  boost::variant<distributed_matrix<double,matrix_row_major>*, distributed_matrix<double,matrix_col_major>*> _ptr;
  wrap_mapping_bc map;
};

} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MATRIX_HPP
