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
    const auto map = solver.default_mapping(dim, rokko::mpi_comm{comm});
    rokko::distributed_crs_matrix mat(map, 3);

    if (map.start_row() == 0) {
      mat.insert(0, {0, 1}, {1., -1.});
    }

    for (int row = std::max(1,map.start_row()); row < std::min(map.end_row(),dim-1); ++row) {
      mat.insert(row, {row-1, row, row+1}, {-1., 2., -1.});
    }

    if (map.end_row() == dim) {
      mat.insert(dim-1, {dim-2, dim-1}, {-1., 2.});
    }

    mat.complete();

    std::vector<std::array<std::string,2>> routines;
    if (library=="anasazi")
      routines = { {"lobpcg", "largest_real"}, {"block_krylov_schur", "largest"}, {"block_davidson", "largest_real"}, {"rtr", "largest_real"} };
    else if (library=="slepc")
      routines = { {"krylovschur", "largest"}, {"lanczos", "largest"}, {"lobpcg", "largest_real"}, {"rqcg", "smallest_real"}, {"subspace", "largest"} };  // excludes not converging "power"

    for (auto routine : routines) {
      std::cout << "routine=" << routine[0] << std::endl;
      auto params_tmp = params;
      params_tmp.set("routine", routine[0]);
      if (routine[0] == "block_davidson")  params_tmp.set("block_size", 10);
      params_tmp.set("wanted_eigenvalues", routine[1]);

      const auto info = solver.diagonalize(mat, params_tmp);

      const auto num_conv = info.get<int>("num_conv");
      if (num_conv == 0)
        throw std::runtime_error("num_conv=0: solver did not converge");

      const auto eigval = solver.eigenvalue(0);
      const auto th_eigval = (routine[0] == "rqcg") ? rokko::laplacian_matrix::eigenvalue(dim, 0)  // smallest one
        : rokko::laplacian_matrix::eigenvalue(dim, dim-1);  // largest one
      EXPECT_NEAR(eigval, th_eigval, eigval*eps);
    }
  } // end if comm != MPI_COMM_NULL
  solver.finalize();
}

TEST(laplacian_crs, eigenvalue) {
  for(auto library : rokko::parallel_sparse_ev::solvers()) {
    run_test(library, MPI_COMM_WORLD);
    run_test(library, create_even_odd_comm_by_split());
    run_test(library, create_even_odd_comm());
    run_test(library, create_even_comm());  // MPI_COMM_NULL for odd rank number
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
