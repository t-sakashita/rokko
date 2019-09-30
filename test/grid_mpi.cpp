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

TEST(grid, grid_processes) {
  MPI_Comm comm = MPI_COMM_WORLD;
  {
    rokko::grid g(comm); // default should be row-major
    // Test public interfaces
    ASSERT_EQ(g.get_comm(), (MPI_Comm)MPI_COMM_WORLD);
    ASSERT_TRUE(g.is_row_major());
    ASSERT_FALSE(g.is_col_major());
  }
  {
    rokko::grid g(comm, rokko::grid_row_major);
    // Test public interfaces
    ASSERT_EQ(g.get_comm(), (MPI_Comm)MPI_COMM_WORLD);
    ASSERT_TRUE(g.is_row_major());
    ASSERT_FALSE(g.is_col_major());

    // Test global values
    if (g.get_nprocs() == 2) {
      ASSERT_EQ(g.get_nprocs(), 2);
      ASSERT_EQ(g.get_nprow(), 1);
      ASSERT_EQ(g.get_npcol(), 2);
    }

    // Test local values
    ASSERT_EQ(g.get_myrow(), g.get_myrank() / g.get_npcol());
    ASSERT_EQ(g.get_mycol(), g.get_myrank() % g.get_npcol());
  }
  {
    rokko::grid g(comm, rokko::grid_col_major);
    // Test public interfaces
    ASSERT_EQ(g.get_comm(), (MPI_Comm)MPI_COMM_WORLD);
    ASSERT_TRUE(g.is_col_major());
    ASSERT_FALSE(g.is_row_major());

    // Test global values
    if (g.get_nprocs() == 2) {
      ASSERT_EQ(g.get_nprocs(), 2);
      ASSERT_EQ(g.get_nprow(), 1);
      ASSERT_EQ(g.get_npcol(), 2);
    }

    // Test local values
    ASSERT_EQ(g.get_myrow(), g.get_myrank() % g.get_nprow());
    ASSERT_EQ(g.get_mycol(), g.get_myrank() / g.get_nprow());
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
