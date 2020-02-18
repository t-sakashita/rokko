/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/laplacian_mfree.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/various_mpi_comm.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-8;

int global_argc;
char** global_argv;

void run_test(MPI_Comm comm) {
  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  if (global_argc >= 2) library_routine = global_argv[1];
  std::string library, routine;
  rokko::split_solver_name(library_routine, library, routine);
  int dim = (global_argc >= 3) ? std::stoi(global_argv[2]) : 100;

  std::cout << "library:routine = " << library_routine << std::endl;

  rokko::parameters params;
  if (!routine.empty()) params.set("routine", routine);
  params.set("block_size", 5);
  params.set("max_iters", 500);
  params.set("conv_tol", eps);
  params.set("wanted_eigenvalues", "largest");

  rokko::parallel_sparse_ev solver(library);
  if (comm != MPI_COMM_NULL) {
    rokko::laplacian_mfree mat(dim, comm);
    rokko::parameters info = solver.diagonalize(mat, params);

    int num_conv = info.get<int>("num_conv");
    if (num_conv == 0)
      throw std::runtime_error("num_conv=0: solver did not converge");

    double eigval = solver.eigenvalue(0);
    double th_eigval = rokko::laplacian_matrix::eigenvalue(dim, dim-1);  // largest one
    EXPECT_NEAR(eigval, th_eigval, eigval*eps);
  }
  solver.finalize();
}

TEST(laplacian_mfree, eigenvalue) {
  run_test(MPI_COMM_WORLD);
  run_test(create_even_odd_comm_by_split());
  run_test(create_even_odd_comm());
  run_test(create_even_comm());
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
