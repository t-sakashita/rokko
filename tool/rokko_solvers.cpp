/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/foreach.hpp>
#include <rokko/solver.hpp>

int main() {
  std::cout << "[serial dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::serial_dense_ev::solvers()) {
    std::cout << "  " << name << std::endl;
  }

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  std::cout << "[parallel dense solvers]\n";
  BOOST_FOREACH(std::string name, rokko::parallel_dense_ev::solvers()) {
    std::cout << "  " << name << std::endl;
  }
#endif // ROKKO_HAVE_PARALLEL_DENSE_SOLVER

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  std::cout << "[parallel sparse solvers]\n";
  BOOST_FOREACH(std::string name, rokko::parallel_sparse_ev::solvers()) {
    std::cout << "  " << name << std::endl;
  }
#endif // ROKKO_HAVE_PARALLEL_SPARSE_SOLVER

  return 0;
}
