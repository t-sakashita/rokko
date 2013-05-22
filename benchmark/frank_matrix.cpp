/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>
#include <ctime>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/config.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  typedef rokko::grid_row_major grid_major;
  //typedef rokko::grid_col_major grid_major;
  //typedef rokko::matrix_col_major matrix_major;
  typedef rokko::matrix_row_major matrix_major;

  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name matrix_size" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  //mkl_set_num_threads(4);

  std::string solver_name(argv[1]);
  unsigned int dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::solver solver(solver_name);
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<grid_major> g(comm);
  int myrank = g.get_myrank();
  int nprocs = g.get_nprocs();

  const int root = 0;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::frank_matrix::generate(mat);

  rokko::localized_vector w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g, solver);

  double start, end;
  try {
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    solver.diagonalize(mat, w, Z);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  std::cout.precision(20);

#ifndef NDEBUG
  if (myrank == root) {
    std::cout << "Computed Eigenvalues= " << w.transpose() << std::endl;
  }
#endif

  if (myrank == 0) {
    double time = end - start;
    //#ifdef _OPENMP_
    std::cout << "num_procs = " << nprocs << std::endl;
    std::cout << "num_threads = " << omp_get_num_threads() << std::endl;
    //std::cout << "num_threads = " << mkl_get_num_threads() << std::endl;
    //#endif
    std::cout << "solver_name = " << solver_name << std::endl;
    std::cout << "matrix = frank" << std::endl;
    std::cout << "dim = " << dim << std::endl;
    std::cout << "time = " << time << std::endl;
    std::cout << "rokko_version = " << ROKKO_VERSION << std::endl;
    std::cout << "hostname = " << boost::asio::ip::host_name() << std::endl;
    std::time_t now = std::time(0);
    std::cout << "date = " << ctime(&now); // << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
