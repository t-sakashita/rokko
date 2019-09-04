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

#ifndef PYROKKO_LOCALIZED_MATRIX_HPP
#define PYROKKO_LOCALIZED_MATRIX_HPP

#include <boost/variant.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/localized_matrix.hpp>
#include <rokko/pyrokko_matrix_major_enum.hpp>

namespace rokko {

namespace py = pybind11;

class wrap_localized_matrix {
public:

  wrap_localized_matrix(matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _ptr = new localized_matrix<double,matrix_col_major>();
    else
      _ptr = new localized_matrix<double,matrix_row_major>();
  }
  
  wrap_localized_matrix(int rows, int cols, matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _ptr = new localized_matrix<double,matrix_col_major>(rows, cols);
    else
      _ptr = new localized_matrix<double,matrix_row_major>(rows, cols);
  }

  localized_matrix<double,matrix_col_major>& col_ver() {
    return *(boost::get<localized_matrix<double,matrix_col_major>*>(_ptr));
  }
  
  localized_matrix<double,matrix_row_major>& row_ver() {
    return *(boost::get<localized_matrix<double,matrix_row_major>*>(_ptr));
  }

  localized_matrix<double,matrix_col_major> const& col_ver() const {
    return *(boost::get<localized_matrix<double,matrix_col_major>*>(_ptr));
  }
  
  localized_matrix<double,matrix_row_major> const& row_ver() const {
    return *(boost::get<localized_matrix<double,matrix_row_major>*>(_ptr));
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

  template <typename MATRIX_MAJOR>
  auto const eigen_ptr() {
    using current_type = localized_matrix<double,MATRIX_MAJOR>;
    return static_cast<typename current_type::super_type*>(boost::get<current_type*>(_ptr));
  }

  py::object get_object() {
    if (is_col)
      return py::cast(eigen_ptr<matrix_col_major>());
    else
      return py::cast(eigen_ptr<matrix_row_major>());
  }

  template <typename MATRIX_MAJOR>
  localized_matrix<double,MATRIX_MAJOR>* get_ptr() {
    return boost::get<localized_matrix<double,MATRIX_MAJOR>*>(_ptr);
  }

  template <typename MATRIX_MAJOR>
  double* get_array_pointer() {
    return get_ptr<MATRIX_MAJOR>()->data();
  }

  void set_ndarray(py::array_t<double> const& mat) {
    py::array_t<double> array;
    if (is_col) {
      array = py::array_t<double>({get_m_local(), get_n_local()}, get_array_pointer<matrix_col_major>(), py::cast(*this));
    } else {
      array = py::array_t<double>({get_m_local(), get_n_local()}, get_array_pointer<matrix_row_major>(), py::cast(*this));
    }

    auto r = array.template mutable_unchecked<2>();
    for (auto i = 0; i < r.shape(0); ++i) {
      for (auto j = 0; j < r.shape(1); ++j) {
        r(i, j) = *mat.data(i, j);
      }
    }
  }

  void set_matrix_col_major(Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> mat) {
    if (!is_col)
      throw std::invalid_argument("Cannot set col-major ndarray to row-major localized_matrix");
    
    col_ver() = static_cast<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(mat);
  }

  void set_matrix_row_major(Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mat) {
    if (is_col)
      throw std::invalid_argument("Cannot set row-major ndarray to col-major localized_matrix");
    
    row_ver() = static_cast<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mat);
  }

private:
  bool is_col;
  boost::variant<localized_matrix<double,matrix_row_major>*, localized_matrix<double,matrix_col_major>*> _ptr;
};

} // end namespace rokko

#endif // PYROKKO_LOCALIZED_MATRIX_HPP
