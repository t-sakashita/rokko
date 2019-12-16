/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/grid.hpp>

#include <gtest/gtest.h>

void run_test(rokko::grid const& g) {
  ASSERT_TRUE(g.get_comm() == MPI_COMM_WORLD);
  // The following tests should be run in single MPI process.
  ASSERT_EQ(g.get_nprocs(), 1);
  ASSERT_EQ(g.get_nprow(), 1);
  ASSERT_EQ(g.get_npcol(), 1);
  ASSERT_EQ(g.get_myrank(), 0);
  ASSERT_EQ(g.get_myrow(), 0);
  ASSERT_EQ(g.get_mycol(), 0);
  ASSERT_EQ(g.calculate_grid_row(g.get_myrank()), 0);
  ASSERT_EQ(g.calculate_grid_col(g.get_myrank()), 0);
}

TEST(grid, grid_1process) {
  MPI_Comm comm = MPI_COMM_WORLD;
  {
    rokko::grid g(comm); // default should be row-major
    ASSERT_TRUE(g.is_row_major());
    ASSERT_FALSE(g.is_col_major());
    run_test(g);
  }
  {
    rokko::grid g(comm, rokko::grid_row_major);
    ASSERT_TRUE(g.is_row_major());
    ASSERT_FALSE(g.is_col_major());
    run_test(g);
  }
  {
    rokko::grid g(comm, rokko::grid_col_major);
    ASSERT_FALSE(g.is_row_major());
    ASSERT_TRUE(g.is_col_major());
    run_test(g);
  }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
