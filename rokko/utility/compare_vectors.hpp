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

template <typename T>
T sign(T a) {
  return (a>0) - (a<0);
}

template <typename MAT1, typename MAT2>
auto norm_diff(const MAT1 source, const MAT2 target) {
  const auto sign_diff = sign(source(0,0) * target(0,0));
  return (source - sign_diff * target).norm();
}

} // namespace rokko

#endif // ROKKO_UTILITY_COMPARE_VECTORS_HPP
