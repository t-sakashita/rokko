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

#include <rokko/mpi_communicator.hpp>

#include <gtest/gtest.h>

TEST(mpi_comm, 2proc) {
  rokko::mpi_comm comm(MPI_COMM_WORLD);
  // Test public interfaces
  ASSERT_TRUE(comm.get_comm() == MPI_COMM_WORLD);

  // Test global values
  ASSERT_EQ(comm.get_nprocs(), 2);
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
