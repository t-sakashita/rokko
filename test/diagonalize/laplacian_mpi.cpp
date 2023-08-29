/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/command_line_parameters.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-5;

int global_argc;
char** global_argv;

template<typename GRID_MAJOR, typename MATRIX_MAJOR>
void test_diagonalize_matrix(std::string const& name) {
  const rokko::mpi_comm comm(MPI_COMM_WORLD);
  constexpr int root_proc = 0;
  if (comm.get_myrank() == root_proc) {
    const std::string grid_major_str = std::is_same_v<GRID_MAJOR,rokko::grid_row_major_t> ? "grid_row_major" : "grid_col_major";
  std::cout << "library=" << name << " grid_major=" << grid_major_str << std::endl;
  }

  constexpr int dim = 10;

  rokko::parallel_dense_ev solver(name);
  if (solver.is_available_grid_major(GRID_MAJOR{})) {
    solver.initialize(global_argc, global_argv);
    const rokko::grid g(comm, GRID_MAJOR{});
    const rokko::mapping_bc<MATRIX_MAJOR> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double, MATRIX_MAJOR> mat(map);
    rokko::laplacian_matrix::generate(mat);
    Eigen::VectorXd w(dim);
    rokko::distributed_matrix<double, MATRIX_MAJOR> Z(map);

    solver.diagonalize(mat, w, Z);

    EXPECT_NEAR(w.array().inverse().sum(), dim * (dim+1) * 0.5, eps);

    solver.finalize();
  }
}

TEST(diagonalize, laplacian_mpi) {
  const auto names = global_argc == 1 ? rokko::parallel_dense_ev::solvers()
    : rokko::get_command_line_args(global_argc, global_argv);

  for(auto const& name : names) {
    test_diagonalize_matrix<rokko::grid_row_major_t,rokko::matrix_col_major>(name);
    test_diagonalize_matrix<rokko::grid_col_major_t,rokko::matrix_col_major>(name);
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
