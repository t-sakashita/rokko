/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/laplacian_mfree.hpp>
#include <rokko/utility/solver_name.hpp>

#include <stdexcept>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  if (argc >= 2) library_routine = argv[1];
  const auto [library, routine] = rokko::split_solver_name(library_routine);
  int dim = (argc >= 3) ? std::stoi(argv[2]) : 100;

  rokko::parameters params;
  if (!routine.empty()) params.set("routine", routine);
  params.set("block_size", 5);
  params.set("max_iters", 500);
  params.set("conv_tol", 1.0e-8);

  rokko::parallel_sparse_ev solver(library);
  rokko::laplacian_mfree mat(dim);
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of Laplacian" << std::endl
	      << "solver = " << library << std::endl
	      << "dimension = " << mat.get_dim() << std::endl;
  
  rokko::parameters info = solver.diagonalize(mat, params);
  
  int num_conv = info.get<int>("num_conv");
  if (num_conv == 0)
    throw std::runtime_error("num_conv=0: solver did not converge");
  std::vector<double> eigvec;
  solver.eigenvector(0, eigvec);
  if (rank == 0) {
    std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
    std::cout << "largest eigenvalues:     ";
    for (int i = 0; i < num_conv; ++i) std::cout << solver.eigenvalue(i) << ' ';
    std::cout << std::endl;
    std::cout << "theoretical eigenvalues: ";
    for (int k = dim - 1; k >= dim - num_conv; --k)
      std::cout << rokko::laplacian_matrix::eigenvalue(dim, k) << ' ';
    std::cout << std::endl;
    std::cout << "largest eigenvector: ";
    for (size_t j = 0; j < eigvec.size(); ++j) std::cout << eigvec[j] << ' ';
    std::cout << std::endl;
  }
  
  solver.finalize();
  MPI_Finalize();
}
