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
#include <rokko/utility/frank_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string solver_name(rokko::parallel_dense_ev::default_solver());
  int dim = 10;
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<int>(argv[2]);

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(5);

  rokko::parallel_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
              << "solver = " << solver_name << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::frank_matrix::generate(mat);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<matrix_major>> mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, 0);

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  //rokko::parameters params;
  //params.set("routine", "aa");
  try {
    solver.diagonalize(mat, eigval, eigvec);
  }
  catch (const char *e) {
    if (myrank == 0) std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<matrix_major>> eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
    std::cout << "largest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
    std::cout << "residual of the largest eigenvalue/vector: |x A x - lambda| = "
              << std::abs(eigvec_loc.col(dim - 1).transpose() * mat_loc * eigvec_loc.col(dim - 1)
                          - eigval(dim - 1))
              << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
