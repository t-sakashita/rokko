/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#define ROKKO_ENABLE_TIMER

#include <mpi.h>
#include <iostream>
#include <fstream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/collective.hpp>

#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>
#include <rokko/config.hpp>
#include <rokko/utility/timer.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  //typedef rokko::matrix_row_major matrix_major;
  typedef rokko::matrix_col_major matrix_major;

  if (argc <= 1) {
    std::cerr << "error: " << argv[0] << " solver_name" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  rokko::timer timer;
  timer.registrate( 1, "diagonalize");
  std::string solver_name(argv[1]);

  int L = 12;
  std::vector<std::pair<int, int> > lattice;
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }
  
  std::cout << "L=" << L << std::endl;

  int dim = 1 << L;
  std::cout << "dim=" << dim << std::endl;
  rokko::parallel_dense_solver solver(solver_name);
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  int myrank = g.get_myrank();
  int nprocs = g.get_nprocs();

  const int root = 0;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  std::cout << "finished generate" << std::endl;

  rokko::localized_vector w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g, solver);

  for (int count=0; count<1; ++count) {
    try {
      solver.diagonalize(mat, w, Z, timer);
    }

    catch (const char *e) {
      std::cout << "Exception : " << e << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 22);
    }
  }
  std::cout << "finished eigensolver" << std::endl;

  /*
  // gather of eigenvectors
  rokko::localized_matrix<matrix_major> eigvec_global;
  rokko::localized_matrix<matrix_major> eigvec_sorted(dim, dim);
  rokko::localized_vector eigval_sorted(dim);
  rokko::gather(Z, eigvec_global, root);
  Z.print();
  //if (myrank == root) {
  //  std::cout << "eigvec:" << std::endl << eigvec_global << std::endl;
  //}

  std::cout.precision(20);
  */

  /*
  std::cout << "w=" << std::endl;
  for (int i=0; i<dim; ++i) {
    std::cout << w[i] << " ";
  }
  std::cout << std::endl;
  */

  if (myrank == 0) {

    std::cout << "num_procs = " << nprocs << std::endl;
#ifdef _OPENMP_
    std::cout << "num_threads = " << omp_get_max_threads() << std::endl;
    //std::cout << "num_threads = " << mkl_get_num_threads() << std::endl;
#endif
    std::cout << "solver_name = " << solver_name << std::endl;
    std::cout << "matrix = frank" << std::endl;
    std::cout << "dim = " << dim << std::endl;
    std::cout << "time = " << timer.get_average(1) << std::endl;
    std::cout << "rokko_version = " << ROKKO_VERSION << std::endl;
    std::cout << "hostname = " << boost::asio::ip::host_name() << std::endl;
    std::time_t now = std::time(0);
    std::cout << "date = " << ctime(&now)<< std::endl;
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
