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
#include <rokko/traits/grid2matrix_major.hpp>
#include <rokko/utility/matrix012.hpp>

#include <gtest/gtest.h>


template<typename GRID_MAJOR>
void run_test(GRID_MAJOR, std::array<int,2> size) {
  constexpr int eigen_major = rokko::eigen3_major<GRID_MAJOR>;
  const rokko::grid g(MPI_COMM_WORLD, size, GRID_MAJOR{});
  ASSERT_TRUE(g.get_comm() == MPI_COMM_WORLD);
  ASSERT_EQ(g.get_nprow() * g.get_npcol(), g.get_nprocs());

  const int lld = (eigen_major == Eigen::RowMajor) ? g.get_npcol() : g.get_nprow();
  ASSERT_EQ(rokko::matrix012::get_index<eigen_major>(lld, g.get_myrow(), g.get_mycol()), g.get_myrank());
}

TEST(grid_matrix, any_shpae) {
  int nprocs, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for (int nprow=1; nprow<=nprocs; ++nprow) {
    if (nprocs % nprow == 0) {
      const int npcol = nprocs / nprow;
      if (myrank == 0)
        std::cout << "nprow=" << nprow << " npcol=" << npcol << std::endl;
      run_test(rokko::grid_row_major, {nprow,npcol});
      run_test(rokko::grid_col_major, {nprow,npcol});
    }
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
