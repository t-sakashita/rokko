/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/helmert_matrix.hpp>
#include <iostream>

using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string library_routine(rokko::parallel_dense_ev::default_solver());
  std::string library, routine;
  int dim = 10;
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = std::stoi(argv[2]);
  rokko::split_solver_name(library_routine, library, routine);

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(std::numeric_limits<unsigned long long>::digits10);

  rokko::parallel_dense_ev solver(library);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of Helmert matrix" << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
      	      << "library:routine = " << library_routine << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_bc<matrix_major> map(dim, 16, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  Eigen::VectorXd diag(dim);
  diag.setLinSpaced(diag.size(), 1, diag.size()); // diag = [1, 2, 3, ..., dim]
  rokko::helmert_matrix::generate_for_given_eigenvalues(mat, diag);
  mat.print();

  Eigen::MatrixXd mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, 0);

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  rokko::parameters params;
  params.set("routine", routine);
  solver.diagonalize(mat, eigval, eigvec, params);

  Eigen::MatrixXd eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
    std::cout << "largest eigenvalues:";
    //for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    for (int i = 0; i < dim; ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
    std::cout << "residual of the largest eigenvalue/vector: |x A x - lambda| = "
              << std::abs(eigvec_loc.col(dim - 1).transpose() * mat_loc * eigvec_loc.col(dim - 1)
                          - eigval(dim - 1))
              << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
