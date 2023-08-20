/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/grid.hpp>
#include <rokko/solver.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/utility/minij_matrix.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template <typename T, typename MATRIX_MAJOR>
void test_trace(rokko::mapping_bc<MATRIX_MAJOR> const& map) {
  constexpr int root_proc = 0;
  auto const& g = map.get_grid();
  const auto dim = map.get_m_global();

  rokko::distributed_matrix<T,MATRIX_MAJOR> mat(map);
  rokko::minij_matrix::generate(mat);
  const auto trace_mat = trace(mat);

  if (g.get_myrank() == root_proc) {
    ASSERT_EQ(trace_mat, static_cast<T>(dim * (dim+1) / 2));
  }
}

TEST(generate_matrix, minij_mpi) {
  constexpr int dim = 1000;
  const rokko::grid g(MPI_COMM_WORLD);

  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    const rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);

    test_trace<double>(map);
    test_trace<float>(map);
    test_trace<int>(map);

    solver.finalize();
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
