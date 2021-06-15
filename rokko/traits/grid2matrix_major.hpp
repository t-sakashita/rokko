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

#pragma once

#include <rokko/eigen3.hpp>
#include <rokko/matrix_major.hpp>
#include <rokko/eigen3/matrix_major.hpp>

namespace rokko {

namespace detail {

template<>
struct eigen3_matrix_major<rokko::grid_row_major_t> {
  static constexpr int value = Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::grid_col_major_t> {
  static constexpr int value = Eigen::ColMajor;
};

} // namespace detail

} // namespace rokko
