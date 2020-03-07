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
#include <rokko/utility/various_mpi_comm.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-8;

int global_argc;
char** global_argv;

void run_test(std::string const& library, MPI_Comm comm) {
  constexpr int dim = 20;
  std::cout << "library = " << library << std::endl;

  rokko::parameters params;
  params.set("num_eigvals", 1);
  params.set("block_size", 5);
  params.set("max_iters", 500);
  params.set("conv_tol", eps);

  rokko::parallel_sparse_ev solver(library);
  if (comm != MPI_COMM_NULL) {
    rokko::laplacian_mfree mat(dim, comm);

    std::vector<std::array<std::string,2>> routines;
    if (library == "anasazi")
      routines = { {"lobpcg", "largest_real"}, {"block_krylov_schur", "largest"}, {"rtr", "largest_real"} };  // excludes not converging "block_davidson"
    else if (library == "slepc")
      routines = { {"krylovschur", "largest"}, {"lanczos", "largest"}, {"lobpcg", "largest_real"}, {"subspace", "largest"} };  // excludes not converging "power"

    for (auto routine : routines) {
      std::cout << "routine=" << routine[0] << std::endl;
      params.set("routine", routine[0]);
      params.set("wanted_eigenvalues", routine[1]);

      rokko::parameters info = solver.diagonalize(mat, params);

      int num_conv = info.get<int>("num_conv");
      if (num_conv == 0)
        throw std::runtime_error("num_conv=0: solver did not converge");

      double eigval = solver.eigenvalue(0);
      double th_eigval = rokko::laplacian_matrix::eigenvalue(dim, dim-1);  // largest one
      EXPECT_NEAR(eigval, th_eigval, eigval*eps);
    }
  } // end if comm != MPI_COMM_NULL
  solver.finalize();
}

TEST(laplacian_mfree, eigenvalue) {
  for(auto library : rokko::parallel_sparse_ev::solvers()) {
    run_test(library, MPI_COMM_WORLD);
    run_test(library, create_even_odd_comm_by_split());
    run_test(library, create_even_odd_comm());
    run_test(library, create_even_comm());  // MPI_COMM_NULL for odd rank number
  }
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
