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

namespace rokko {

namespace detail {
    
template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct eigen3_matrix_major;

template<>
struct eigen3_matrix_major<rokko::matrix_row_major> {
  static const int value= Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major> {
  static const int value = Eigen::ColMajor;
};

} // end namespace detail
  
template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct localized_matrix : public Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, detail::eigen3_matrix_major<MATRIX_MAJOR>::value > {
public:
  localized_matrix() : Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, detail::eigen3_matrix_major<MATRIX_MAJOR>::value >() {};
  localized_matrix(int rows, int cols) : Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, detail::eigen3_matrix_major<MATRIX_MAJOR>::value >(rows, cols) {};
};

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
