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

namespace rokko {

// Assume that n >= 1
int find_power_of_two(int n) {
  int p = 0;
  while (n > 1) {
    n /= 2;
    ++p;
  };
  return p;
}

} // namespace rokko
