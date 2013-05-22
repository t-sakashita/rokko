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

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/frank_matrix.hpp>

#define ROKKO_ENABLE_TIMER
#define ROKKO_ENABLE_TIMER_TRACE

#include <rokko/utility/timer.hpp>

int main(int argc, char **argv) {
  int dim;
  if (argc > 1) {
    dim = boost::lexical_cast<int>(argv[1]);
  } else {
    std::cin >> dim;
  }

  rokko::timer timer;
  timer.registrate( 1, "generate");
  timer.registrate( 2, "eigenvalue");

  timer.start(1);
  rokko::localized_matrix<rokko::matrix_col_major> mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "dimension = " << dim << std::endl;
  std::cout << "[elements of frank matrix]" << std::endl;
  std::cout << mat << std::endl;
  timer.stop(1);

  timer.start(2);
  std::cout << "[eigenvalues of frank matrix]" << std::endl;
  double sum = 0;
  for (int i = 0; i < dim; ++i) {
    double ev = rokko::frank_matrix::eigenvalue(dim, i);
    sum += ev;
    std::cout << ev << "  ";
  }
  std::cout << std::endl;
  timer.stop(2);

  std::cout << "[sum of eigenvalues of frank matrix]" << std::endl;
  std::cout << sum << std::endl;
  
  timer.summarize();
  // std::cout << timer.get_average(1) << std::endl;
}
