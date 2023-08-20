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
#include <rokko/eigen3.hpp>
#include <rokko/traits/grid2matrix_major.hpp>
#include <rokko/utility/matrix012.hpp>

#include <gtest/gtest.h>


template<typename GRID_MAJOR>
void run_test(GRID_MAJOR) {
  const rokko::grid g(MPI_COMM_WORLD, GRID_MAJOR{});
  ASSERT_TRUE(g.get_comm() == MPI_COMM_WORLD);
  ASSERT_EQ(g.get_nprow() * g.get_npcol(), g.get_nprocs());

  Eigen::Matrix<int, Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<GRID_MAJOR>> mat(g.get_nprow(), g.get_npcol());
  rokko::matrix012::generate(mat);
  ASSERT_EQ(mat(g.get_myrow(), g.get_mycol()), g.get_myrank());
}

TEST(grid_matrix, default_shape) {
  // default shape: square or nearly square
  run_test(rokko::grid_row_major);
  run_test(rokko::grid_col_major);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
