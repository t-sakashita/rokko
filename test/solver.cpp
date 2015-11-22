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

#define BOOST_TEST_MODULE test_solver
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_solver) {
  BOOST_FOREACH(std::string name, rokko::serial_dense_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    solver.finalize();
  }

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  BOOST_FOREACH(std::string name, rokko::parallel_dense_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    solver.finalize();
  }
#endif // ROKKO_HAVE_PARALLEL_DENSE_SOLVER

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  BOOST_FOREACH(std::string name, rokko::parallel_sparse_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::parallel_sparse_ev solver(name);
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    solver.finalize();
  }
#endif // ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
}
