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

#include <mpi.h>
#include <iostream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/collective.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

#include <boost/lexical_cast.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  //typedef rokko::matrix_row_major matrix_major;
  typedef rokko::matrix_col_major matrix_major;

  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name matrix_size" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  std::string solver_name(argv[1]);
  unsigned int dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::parallel_dense_solver solver(solver_name);
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  int myrank = g.get_myrank();

  const int root = 1;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::frank_matrix::generate(mat);
  mat.print();

  rokko::localized_matrix<matrix_major> lmat(dim, dim);
  rokko::gather(mat, lmat, root);
  if (myrank == root)
    std::cout << "lmat:" << std::endl << lmat << std::endl;

  solver.finalize();
  MPI_Finalize();
  return 0;
}
