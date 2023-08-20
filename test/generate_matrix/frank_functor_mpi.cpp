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

class frank {
public:
  frank(int dim_in) :  dim(dim_in) {}
  double operator()(int i, int j) {
    return dim - std::max(i, j);
  }
private:
  int dim;
};

TEST(distributed_matrix, frank_functor_mpi) {
  constexpr int dim = 10;
  rokko::grid g(MPI_COMM_WORLD);

  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
    mat.generate(frank(dim));

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
