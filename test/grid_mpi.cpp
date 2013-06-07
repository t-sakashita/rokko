/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

#include <rokko/grid.hpp>

#define BOOST_TEST_MODULE test_grid_mpi
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_grid_mpi) {
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    rokko::grid g(comm, rokko::grid_col_major);
    // Test public interfaces
    BOOST_CHECK_EQUAL(g.get_comm(), MPI_COMM_WORLD);
    BOOST_CHECK( g.is_col_major());
    BOOST_CHECK(!g.is_row_major());

    // Test global values
    BOOST_CHECK_EQUAL(g.get_nprocs(), 2);
    BOOST_CHECK_EQUAL(g.get_nprow(), 1);
    BOOST_CHECK_EQUAL(g.get_npcol(), 2);

    // Test local values
    BOOST_CHECK_EQUAL(g.get_myrow(), 0);
    BOOST_CHECK_EQUAL(g.get_mycol(), g.get_myrank());
    MPI_Finalize();
}
