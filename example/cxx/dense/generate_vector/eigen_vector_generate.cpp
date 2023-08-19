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

#include <rokko/eigen3/generate_vector.hpp>
#include <iostream>

int main(int argc, char *argv[]) {
  const auto f = [](int i) { return static_cast<double>(2*i); };

  constexpr unsigned int dim = 10;
  Eigen::VectorXd vec(dim);
  rokko::generate(vec, f);

  std::cout << "vec:\n"
            << vec.transpose() << std::endl;
}
