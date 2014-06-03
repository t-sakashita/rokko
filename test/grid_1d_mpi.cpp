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

#include <rokko/grid_1d.hpp>

#define BOOST_TEST_MODULE test_grid_1d_mpi
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_grid_1d_mpi) {
  MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
           &boost::unit_test::framework::master_test_suite().argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid_1d g(comm);
  // Test public interfaces
  BOOST_CHECK_EQUAL(g.get_comm(), (MPI_Comm)MPI_COMM_WORLD);

  // Test global values
  BOOST_CHECK_EQUAL(g.get_nprocs(), 2);
  MPI_Finalize();
}
