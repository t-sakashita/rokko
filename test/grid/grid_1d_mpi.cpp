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

#include <rokko/grid_1d.hpp>

#include <gtest/gtest.h>

TEST(grid, grid_1d_mpi) {
  const MPI_Comm comm = MPI_COMM_WORLD;
  const rokko::grid_1d g(comm);
  // Test public interfaces
  ASSERT_TRUE(g.get_comm() == MPI_COMM_WORLD);

  // Test global values
  ASSERT_EQ(g.get_nprocs(), 2);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
