/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/foreach.hpp>
#include <rokko/solver_factory.hpp>

int main() {
  std::cout << "[serial dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::solver_factory::serial_dense_solver_names()) {
    std::cout << "  " << name << std::endl;
  }
  std::cout << "[parallel dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::solver_factory::parallel_dense_solver_names()) {
    std::cout << "  " << name << std::endl;
  }
  return 0;
}
