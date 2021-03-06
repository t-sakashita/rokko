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
#include <rokko/utility/solver_name.hpp>

#include <stdexcept>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  std::string library, routine;
  if (argc >= 2) library_routine = argv[1];
  rokko::split_solver_name(library_routine, library, routine);
  
  int dim = (argc >= 3) ? std::stoi(argv[2]) : 100;

  rokko::parameters params;
  if (!routine.empty()) params.set("routine", routine);
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  rokko::parallel_sparse_ev solver(library);
  auto map = solver.default_mapping(dim, rokko::mpi_comm{MPI_COMM_WORLD});
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
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of Laplacian" << std::endl
	      << "solver = " << library << std::endl
	      << "dimension = " << dim << std::endl;
  
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
    for (int j = 0; j < eigvec.size(); ++j) std::cout << eigvec[j] << ' ';
    std::cout << std::endl;
  }
  
  solver.finalize();
  MPI_Finalize();
}
