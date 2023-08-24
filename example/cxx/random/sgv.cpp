/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <random>

using vector_t = Eigen::VectorXd;

int main() {
  constexpr int n = 6;
  vector_t v(n);
  std::mt19937 engine(12345lu);
  std::normal_distribution<> dist(0.0, 1.0);
  for (int i = 0; i < n; ++i) v(i) = dist(engine);
  std::cout << "v: " << v << std::endl;
}
