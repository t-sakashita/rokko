/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_COMPARE_VECTORS_HPP
#define ROKKO_UTILITY_COMPARE_VECTORS_HPP

#include <rokko/eigen3.hpp>

namespace rokko {

template<typename T, int SIZE>
bool is_equal_signed(Eigen::Ref<Eigen::Vector<T, SIZE>> source, Eigen::Ref<Eigen::Vector<T, SIZE>> target) {
  int sign;
  if (source(0) - target(0) < EPSS_QR)  sign = 1;
  else if (source(0) + target(0) < EPSS_QR)  sign = -1;
  else return false;
  
  return source.all() == sign*target.all();
}

} // namespace rokko

#endif // ROKKO_UTILITY_COMPARE_VECTORS_HPP
