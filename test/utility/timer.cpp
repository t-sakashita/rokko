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

#include <iostream>
#include <rokko/eigen3.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/timer.hpp>

int main(int argc, char **argv) {
  int dim;
  if (argc > 1) {
    dim = std::stoi(argv[1]);
  } else {
    dim = 1000;
  }

  rokko::timer timer;
  timer.registrate(1, "generate");
  timer.registrate(2, "eigenvalue");

  rokko::global_timer::registrate(1, "generate");
  rokko::global_timer::registrate(2, "eigenvalue");
  timer.start(1);
  rokko::global_timer::start(1);
  Eigen::MatrixXd mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "dimension = " << dim << std::endl;
  std::cout << "[elements of frank matrix]" << std::endl;
  std::cout << mat << std::endl;
  timer.stop(1);
  rokko::global_timer::stop(1);

  timer.start(2);
  rokko::global_timer::start(2);
  std::cout << "[eigenvalues of frank matrix]" << std::endl;
  double sum = 0;
  for (int i = 0; i < dim; ++i) {
    double ev = rokko::frank_matrix::eigenvalue(dim, i);
    sum += ev;
    std::cout << ev << "  ";
  }
  std::cout << std::endl;
  timer.stop(2);
  rokko::global_timer::stop(2);

  std::cout << "[sum of eigenvalues of frank matrix]" << std::endl;
  std::cout << sum << std::endl;
  
  timer.summarize();
  rokko::global_timer::summarize();
}
