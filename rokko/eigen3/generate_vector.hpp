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

#pragma once

#include <iostream>
#include <rokko/eigen3.hpp>

namespace rokko {

template<typename T, int ROWS, class FUNC>
void generate(Eigen::Vector<T,ROWS>& vec, FUNC func) {
  for(int i = 0; i < vec.size(); ++i) {
    vec(i) = func(i);
  }
}

template<typename T, int ROWS>
void generate(Eigen::Vector<T,ROWS>& vec, std::function<T(int)> const& func) {
  for(int i = 0; i < vec.size(); ++i) {
    vec(i) = func(i);
  }
}

} // namespace rokko
