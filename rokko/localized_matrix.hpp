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

#include <iostream>
#include <rokko/eigen3.hpp>

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

template<typename MATRIX_MAJOR>
constexpr int eigen3_major = detail::eigen3_matrix_major<MATRIX_MAJOR>::value;

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
};


template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, FUNC func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, class FUNC>
void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, std::function<T(int, int)> const& func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
