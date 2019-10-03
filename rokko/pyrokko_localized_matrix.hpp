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
      _ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>();
    else
      _ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>();
  }
  
  wrap_localized_matrix(int rows, int cols, matrix_major_enum const& major = col) : is_col(major == col) {
    if (is_col)
      _ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>(rows, cols);
    else
      _ptr = new Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(rows, cols);
  }

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& col_ver() {
    return *(boost::get<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(_ptr));
  }
  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& row_ver() {
    return *(boost::get<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(_ptr));
  }

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> const& col_ver() const {
    return *(boost::get<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(_ptr));
  }
  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> const& row_ver() const {
    return *(boost::get<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(_ptr));
  }

  void generate(std::function<double(int, int)> const& func) {
    if (is_col)
      rokko::generate(col_ver(), func);
    else
      rokko::generate(row_ver(), func);
  }

  bool is_major_col() const {
    return is_col;
  }
  
  std::string get_major_string() const {
    return is_col ? "col" : "row";
  }

  template <int MATRIX_MAJOR>
  auto const eigen_ptr() {
    using current_type = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>;
    return boost::get<current_type*>(_ptr);
  }

  py::object get_object() {
    if (is_col)
      return py::cast(eigen_ptr<Eigen::ColMajor>());
    else
      return py::cast(eigen_ptr<Eigen::RowMajor>());
  }

  template <int MATRIX_MAJOR>
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>* get_ptr() {
    return boost::get<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>*>(_ptr);
  }

  template <int MATRIX_MAJOR>
  double* get_array_pointer() {
    return get_ptr<MATRIX_MAJOR>()->data();
  }

private:
  bool is_col;
  boost::variant<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*> _ptr;
};

} // end namespace rokko

#endif // PYROKKO_LOCALIZED_MATRIX_HPP
