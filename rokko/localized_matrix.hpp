/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LOCALIZED_MATRIX_HPP
#define ROKKO_LOCALIZED_MATRIX_HPP

#include <rokko/matrix_major.hpp>
#include <Eigen/Dense>
#include <iostream>

namespace rokko {

namespace detail {
    
template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct eigen3_matrix_major;

template<>
struct eigen3_matrix_major<rokko::matrix_row_major> {
  static const int value = Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major> {
  static const int value = Eigen::ColMajor;
};

} // end namespace detail
  
template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct localized_matrix : public Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, detail::eigen3_matrix_major<MATRIX_MAJOR>::value > {
public:
  typedef MATRIX_MAJOR major_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, detail::eigen3_matrix_major<major_type>::value>
          super_type;
  typedef localized_matrix<major_type> matrix_type;

  localized_matrix() : super_type() {};
  localized_matrix(int rows, int cols) : super_type(rows, cols) {};

  template <typename T>
  localized_matrix(T const& other) : super_type(other) {}; 
  template <typename T>
  matrix_type& operator=(T const& other) { super_type::operator=(other); return *this; } 

  int translate_l2g_row(const int& local_i) const {
    return local_i;
  }

  int translate_l2g_col(const int& local_j) const {
    return local_j;
  }

  int translate_g2l_row(const int& global_i) const {
    return global_i;
  }

  int translate_g2l_col(const int& global_j) const {
    return global_j;
  }

  int get_m_global() const { return this->rows(); }
  int get_n_global() const { return this->cols(); }

  int get_m_local() const { return this->rows(); }
  int get_n_local() const { return this->cols(); }

  bool is_gindex_myrow(const int& global_i) const {
    return true;
  }

  bool is_gindex_mycol(const int& global_j) const {
    return true;
  }

  bool is_gindex(const int& global_i, const int& global_j) const {
    return true;
  }

  void set_local(int local_i, int local_j, double value) {
    this->operator()(local_i, local_j) = value;
  }

  double get_local(int local_i, int local_j) const {
    return this->operator()(local_i, local_j);
  }
  
  void update_local(int local_i, int local_j, double value) {
    this->operator()(local_i, local_j) += value;
  }

  void set_zeros() {
    this->setZero();
  }

  bool is_row_major() const {
    return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value;
  }
  bool is_col_major() const {
    return boost::is_same<MATRIX_MAJOR, matrix_col_major>::value;
  }

  void print() const {
    std::cout << *this << std::endl;
  }

};

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
