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
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

int dim_global;

double frank_calculate_matrix_element(int i, int j) {
  return dim_global - std::max(i, j);
}

template<typename GRID_MAJOR, typename MATRIX_MAJOR>
void test_generate_matrix(std::string const& name) {
  const rokko::mpi_comm comm(MPI_COMM_WORLD);
  constexpr int root_proc = 0;
  if (comm.get_myrank() == root_proc) {
    const std::string grid_major_str = std::is_same_v<GRID_MAJOR,rokko::grid_row_major_t> ? "grid_row_major" : "grid_col_major";
    std::cout << "library=" << name << " grid_major=" << grid_major_str << std::endl;
  }

  constexpr int dim = 10;
  dim_global = dim;

  rokko::parallel_dense_ev solver(name);
  if (solver.is_available_grid_major(GRID_MAJOR{})) {
    solver.initialize(global_argc, global_argv);
    const rokko::grid g(comm, GRID_MAJOR{});
    const rokko::mapping_bc<MATRIX_MAJOR> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double,MATRIX_MAJOR> mat(map);
    mat.generate(&frank_calculate_matrix_element);

    const auto dim_proc = (g.get_myrank() == root_proc) ? dim : 0;
    Eigen::MatrixXd lmat_gather(dim_proc, dim_proc);
    rokko::gather(mat, lmat_gather, root_proc);

    if (comm.get_myrank() == root_proc) {
      Eigen::MatrixXd lmat(dim, dim);
      rokko::frank_matrix::generate(lmat);
      ASSERT_TRUE(lmat_gather == lmat);
    }

    mat.print();
    solver.finalize();
  }
}

TEST(distributed_matrix, frank_functor_mpi) {
  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    test_generate_matrix<rokko::grid_row_major_t,rokko::matrix_col_major>(name);
    test_generate_matrix<rokko::grid_col_major_t,rokko::matrix_col_major>(name);
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
