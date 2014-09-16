/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>,
*                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/foreach.hpp>
#include <rokko/solver.hpp>
#include <rokko/parallel_sparse_solver.hpp>

int main() {
  std::cout << "[serial dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::serial_dense_solver::solvers()) {
    std::cout << "  " << name << std::endl;
  }
  std::cout << "[parallel dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::parallel_dense_solver::solvers()) {
    std::cout << "  " << name << std::endl;
  }
  //#if defined(BUILD_ANASAZI) || defined(BUILD_SLEPC)
  std::cout << "[parallel sparse solvers]\n";
  BOOST_FOREACH(std::string name, rokko::parallel_sparse_solver::solvers()) {
    std::cout << "  " << name << std::endl;
  }
  //#endif // defined(BUILD_ANASAZI) || defined(BUILD_SLEPC)
  return 0;
}
