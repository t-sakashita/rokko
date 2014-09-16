/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>,
*               2014-2014 by Synge Todo <wistaria@comp-phys.org>,
*               2014-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

//#include <rokko/solver.hpp>
#include <rokko/parallel_sparse_solver.hpp>

#define BOOST_TEST_MODULE test_parallel_sparse_solvers
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_parallel_sparse_solvers) {
  BOOST_FOREACH(std::string name, rokko::parallel_sparse_solver::solvers()) {
    std::cerr << name << std::endl;
  }
}
