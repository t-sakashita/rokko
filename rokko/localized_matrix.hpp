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

#ifndef ROKKO_LOCALIZED_MATRIX_HPP
#define ROKKO_LOCALIZED_MATRIX_HPP

#include <rokko/eigen3.hpp>
#include <iostream>
#include <boost/type_traits.hpp>

namespace rokko {

namespace detail {
    
template<typename MATRIX_MAJOR>
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
  
template<typename T, typename MATRIX_MAJOR = rokko::matrix_col_major,
         int ROWS = Eigen::Dynamic, int COLS = Eigen::Dynamic>
class localized_matrix : public Eigen::Matrix<T, ROWS, COLS,
  detail::eigen3_matrix_major<MATRIX_MAJOR>::value> {
public:
  typedef T value_type;
  typedef MATRIX_MAJOR major_type;
  typedef Eigen::Matrix<value_type, ROWS, COLS,
    detail::eigen3_matrix_major<major_type>::value> super_type;
  typedef localized_matrix<value_type, major_type, ROWS, COLS> matrix_type;

  localized_matrix() : super_type() {}
  localized_matrix(int rows, int cols) : super_type(rows, cols) {}

  template<typename U>
  localized_matrix(U const& other) : super_type(other) {}
  template<typename U>
  matrix_type& operator=(U const& other) { super_type::operator=(other); return *this; } 

  int translate_l2g_row(int local_i) const { return local_i; }
  int translate_l2g_col(int local_j) const { return local_j; }
  int translate_g2l_row(int global_i) const { return global_i; }
  int translate_g2l_col(int global_j) const { return global_j; }

  int get_m_global() const { return super_type::rows(); }
  int get_n_global() const { return super_type::cols(); }

  int get_m_local() const { return super_type::rows(); }
  int get_n_local() const { return super_type::cols(); }

  bool is_gindex_myrow(int) const { return true; }
  bool is_gindex_mycol(int) const { return true; }
  bool is_gindex(int, int) const { return true; }

  void set_local(int local_i, int local_j, value_type value) {
    this->operator()(local_i, local_j) = value;
  }
  void update_local(int local_i, int local_j, value_type value) {
    this->operator()(local_i, local_j) += value;
  }
  value_type get_local(int local_i, int local_j) const {
    return this->operator()(local_i, local_j);
  }
  
  void set_global(int global_i, int global_j, value_type value) {
    set_local(global_i, global_j, value);
  }
  void update_global(int global_i, int global_j, value_type value) {
    update_local(global_i, global_j, value);
  }
  value_type get_global(int global_i, int global_j) {
    return get_local(global_i, global_j);
  }
  value_type get_global_checked(int global_i, int global_j) {
    return get_local(global_i, global_j);
  }

  template<class FUNC>
  void generate(FUNC func) {
    for(int local_i = 0; local_i < get_m_local(); ++local_i) {
      for(int local_j = 0; local_j < get_n_local(); ++local_j) {
        set_local(local_i, local_j, func(local_i, local_j));
      }
    }
  }

  void set_zeros() { super_type::setZero(); }

  bool is_row_major() const { return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value; }
  bool is_col_major() const { return boost::is_same<MATRIX_MAJOR, matrix_col_major>::value; }

  void print() const { std::cout << *this << std::endl; }
};

typedef localized_matrix<float> flmatrix;
typedef localized_matrix<double> dlmatrix;
typedef localized_matrix<std::complex<float> > clmatrix;
typedef localized_matrix<std::complex<double> > zlmatrix;

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
