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

#include <rokko/grid.hpp>
#include <rokko/solver.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template<typename T, typename MATRIX_MAJOR>
void eigen_2_distributed(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>>& lmat, rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
  auto const& map = mat.get_mapping();
  if (map.get_m_global() != map.get_n_global())
    throw std::invalid_argument("frank_matrix::generate() : non-square matrix");
  for(int local_i = 0; local_i < map.get_m_local(); ++local_i) {
    for(int local_j = 0; local_j < map.get_n_local(); ++local_j) {
      const auto global_i = map.translate_l2g_row(local_i);
      const auto global_j = map.translate_l2g_col(local_j);
      mat.set_local(local_i, local_j, lmat(global_i, global_j));
    }
  }
}

TEST(eigen2distributed_matrix, eigen2distributed_matrix) {
  constexpr int dim = 10;
  const rokko::grid g(MPI_COMM_WORLD);

  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    const rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
    Eigen::MatrixXd lmat(dim, dim);
    rokko::frank_matrix::generate(lmat);
    eigen_2_distributed(lmat, mat);
    rokko::frank_matrix::generate(mat);

    constexpr int root_proc = 0;
    const auto dim_proc = (g.get_myrank() == root_proc) ? dim : 0;
    Eigen::MatrixXd lmat_gather(dim_proc, dim_proc);
    rokko::gather(mat, lmat_gather, root_proc);

    if (g.get_myrank() == root_proc) {
      Eigen::MatrixXd lmat(dim, dim);
      rokko::frank_matrix::generate(lmat);
      ASSERT_TRUE(lmat_gather == lmat);
    }

    mat.print();
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
