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

constexpr int root_proc = 0;

template <typename T, typename MATRIX_MAJOR>
void test_trace(rokko::mapping_bc<MATRIX_MAJOR> const& map) {
  rokko::distributed_matrix<T,MATRIX_MAJOR> mat(map);
  rokko::minij_matrix::generate(mat);
  const auto trace_mat = trace(mat);

  auto const& g = map.get_grid();
  if (g.get_myrank() == root_proc) {
    const auto dim = map.get_m_global();
    ASSERT_EQ(trace_mat, static_cast<T>(dim * (dim+1) / 2));
  }
}

template<typename GRID_MAJOR, typename MATRIX_MAJOR>
void run_test(std::string const& name) {
  const rokko::mpi_comm comm(MPI_COMM_WORLD);
  constexpr int root_proc = 0;
  if (comm.get_myrank() == root_proc) {
    const std::string grid_major_str = std::is_same_v<GRID_MAJOR,rokko::grid_row_major_t> ? "grid_row_major" : "grid_col_major";
    std::cout << "library=" << name << " grid_major=" << grid_major_str << std::endl;
  }

  constexpr int dim = 1000;

  rokko::parallel_dense_ev solver(name);
  if (solver.is_available_grid_major(GRID_MAJOR{})) {
    solver.initialize(global_argc, global_argv);
    const rokko::grid g(MPI_COMM_WORLD, GRID_MAJOR{});
    const rokko::mapping_bc<MATRIX_MAJOR> map = solver.default_mapping(dim, g);

    test_trace<double>(map);
    test_trace<float>(map);
    test_trace<int>(map);

    solver.finalize();
  }
}

TEST(generate_matrix, minij_mpi) {
  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    run_test<rokko::grid_row_major_t,rokko::matrix_col_major>(name);
    run_test<rokko::grid_col_major_t,rokko::matrix_col_major>(name);
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
