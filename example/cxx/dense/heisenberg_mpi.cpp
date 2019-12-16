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
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string solver_name(rokko::parallel_dense_ev::default_solver());
  int L = 8;
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) L = boost::lexical_cast<int>(argv[2]);

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(5);

  int dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1) % L));
  }

  rokko::parallel_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
              << "solver = " << solver_name << std::endl
              << "L = " << L << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  Eigen::MatrixXd mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, 0);

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  try {
    solver.diagonalize(mat, eigval, eigvec);
  } catch (const char *e) {
    if (myrank == 0) std::cerr << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  Eigen::MatrixXd eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
    std::cout << std::endl;
    std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
              << std::abs(eigvec_loc.col(0).transpose() * mat_loc * eigvec_loc.col(0) - eigval(0))
              << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
