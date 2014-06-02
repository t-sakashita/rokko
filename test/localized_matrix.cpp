/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Yuichi Motoyama <y-motoyama@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#define BOOST_TEST_MODULE test_distributed_matrix
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_localized_matrix) {

  // M = 1 2 3
  //     4 5 6
  //     7 8 9
  //
  // a = 5
  // 
  // u = 1
  //     2
  //     3
  //
  // v = 4
  //     5
  //     6
  // 
  // w = au + Mv = 37
  //               87
  //              137

  int dim = 3;
  rokko::localized_matrix<> M(dim,dim);
  M << 1,2,3,4,5,6,7,8,9;
  double a = 5.0;
  rokko::localized_vector  u(dim);
  u << 1,2,3;
  rokko::localized_vector v(dim);
  v << 4,5,6;

  rokko::localized_vector w = a*u+M*v;
  BOOST_CHECK_EQUAL(w[0],37.);
  BOOST_CHECK_EQUAL(w[1],87.);
  BOOST_CHECK_EQUAL(w[2],137.);

  // u^t M = 30 36 42
  BOOST_CHECK_EQUAL( (u.transpose()*M.col(0)), 30.);
  BOOST_CHECK_EQUAL( (u.transpose()*M.col(1)), 36.);
  BOOST_CHECK_EQUAL( (u.transpose()*M.col(2)), 42.);
}
