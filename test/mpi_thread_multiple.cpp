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

#define BOOST_TEST_MODULE test_mpi_thread_multiple
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include <mpi.h>
#include <iostream>

BOOST_AUTO_TEST_CASE(test_mpi_thread_multiple) {
  const int required = MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&boost::unit_test::framework::master_test_suite().argc,
                  &boost::unit_test::framework::master_test_suite().argv,
                  required, &provided);
  BOOST_CHECK_EQUAL(provided, required);
  if (provided == required) {
    std::cerr << "MPI_THREAD_MULTIPLE is supported\n";
  } else {
    std::cerr << "MPI_THREAD_MULTIPLE is NOT supported\n";
  }
  MPI_Finalize();
}
