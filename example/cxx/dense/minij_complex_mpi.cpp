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
#include <rokko/mapping_bc.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/minij_matrix.hpp>
#include <iostream>


using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  const MPI_Comm comm = MPI_COMM_WORLD;
  const std::string library_routine = (argc >= 2) ? argv[1] : rokko::parallel_dense_ev::default_solver();
  const int dim = (argc >= 3) ? std::stoi(argv[2]) : 10;
  const auto [library, routine] = rokko::split_solver_name(library_routine);

  const rokko::grid g(comm);
  const auto myrank = g.get_myrank();

  std::cout.precision(5);

  rokko::parallel_dense_ev solver(library);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of minij matrix" << std::endl
	      << "library:routine = " << library_routine << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
	      << "routine = " << routine << std::endl
              << "dimension = " << dim << std::endl;

  const rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<std::complex<double>, matrix_major> mat(map);
  rokko::minij_matrix::generate(mat);
  Eigen::MatrixXcd mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, 0);

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<std::complex<double>, matrix_major> eigvec(map);
  rokko::parameters params;
  params.set("routine", routine);
  solver.diagonalize(mat, eigval, eigvec, params);

  Eigen::MatrixXcd eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
    bool sorted = true;
    for (int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
    if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";

    std::cout << "largest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
    std::cout << "eigenvectors:\n";
    std::cout << eigvec_loc << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
