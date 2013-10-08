/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/frank_matrix.hpp>

#define BOOST_TEST_MODULE test_scatter
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include <boost/random.hpp>

template<typename GRID_MAJOR, typename DIST_MAT_MAJOR, typename LOC_MAT_MAJOR>
bool run_test(MPI_Comm comm, int dim, GRID_MAJOR const& grid_major, DIST_MAT_MAJOR const&, LOC_MAT_MAJOR const&) {
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // seme seed for all processes
  boost::mt19937 eng;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    rng(eng, boost::uniform_real<>());
  
  rokko::localized_matrix<LOC_MAT_MAJOR> lmat(dim, dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      lmat(i, j) = rng();
    }
  }
#ifndef NDEBUG
  if (rank == 0) std::cout << lmat << std::endl;
#endif
  
  int success_local = 1;
  rokko::solver solver(rokko::solver_factory::instance()->solver_names()[0]);
  rokko::grid g(comm, grid_major);
  for (int r = 0; r < size; ++r) {
    rokko::distributed_matrix<DIST_MAT_MAJOR> mat(dim, dim, g, solver);
    rokko::scatter(lmat, mat, r);
#ifndef NDEBUG
    mat.print();
#endif
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        if (mat.is_gindex(i, j) && std::abs(mat.get_global(i, j) - lmat(i, j)) >  10e-12)
          success_local = 0;
      }
    }
  }

  int success;
  MPI_Allreduce(&success_local, &success, 1, MPI_INT, MPI_PROD, comm);
  return success;
}

BOOST_AUTO_TEST_CASE(test_scatter) {
  MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (boost::unit_test::framework::master_test_suite().argc > 1) {
    dim = boost::lexical_cast<int>(boost::unit_test::framework::master_test_suite().argv[1]);
  }
  if (rank == 0) std::cout << "dimension = " << dim << std::endl;

  if (rank == 0) std::cout << "test for grid = col, dist = col, loc = col\n";
  bool success_col_col_col =
    run_test(comm, dim, rokko::grid_col_major, rokko::matrix_col_major(), rokko::matrix_col_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_col_col_col, true);

  /*
  if (rank == 0) std::cout << "test for grid = col, dist = col, loc = row\n";
  bool success_col_col_row =
    run_test(comm, dim, rokko::grid_col_major, rokko::matrix_col_major(), rokko::matrix_row_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_col_col_row, true);

  if (rank == 0) std::cout << "test for grid = col, dist = row, loc = col\n";
  bool success_col_row_col =
    run_test(comm, dim, rokko::grid_col_major, rokko::matrix_row_major(), rokko::matrix_col_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_col_row_col, true);

  if (rank == 0) std::cout << "test for grid = col, dist = row, loc = row\n";
  bool success_col_row_row =
    run_test(comm, dim, rokko::grid_col_major, rokko::matrix_row_major(), rokko::matrix_row_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_col_row_row, true);
  */

  if (rank == 0) std::cout << "test for grid = row, dist = col, loc = col\n";
  bool success_row_col_col =
    run_test(comm, dim, rokko::grid_row_major, rokko::matrix_col_major(), rokko::matrix_col_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_row_col_col, true);

  /*
  if (rank == 0) std::cout << "test for grid = row, dist = col, loc = row\n";
  bool success_row_col_row =
    run_test(comm, dim, rokko::grid_row_major, rokko::matrix_col_major(), rokko::matrix_row_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_row_col_row, true);

  if (rank == 0) std::cout << "test for grid = row, dist = row, loc = col\n";
  bool success_row_row_col =
    run_test(comm, dim, rokko::grid_row_major, rokko::matrix_row_major(), rokko::matrix_col_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_row_row_col, true);

  if (rank == 0) std::cout << "test for grid = row, dist = row, loc = row\n";
  bool success_row_row_row =
    run_test(comm, dim, rokko::grid_row_major, rokko::matrix_row_major(), rokko::matrix_row_major());
  if (rank == 0) BOOST_CHECK_EQUAL(success_row_row_row, true);
  */

  MPI_Finalize();
}
