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
#include <rokko/rokko.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <boost/lexical_cast.hpp>

typedef rokko::matrix_col_major matrix_major;
// typedef rokko::matrix_row_major matrix_major;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  unsigned int dim = 10;
  std::string solver_name(rokko::parallel_dense_solver::default_solver());
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::grid g;
  if (g.get_myrank() == 0) std::cout << "dimension = " << dim << std::endl;
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());

  rokko::parallel_dense_solver solver(solver_name);
  solver.initialize(argc, argv);
  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::frank_matrix::generate(mat);
  mat.print();

  rokko::localized_matrix<matrix_major> lmat(dim, dim);
  rokko::gather(mat, lmat, 0);
  if (g.get_myrank() == 0)
    std::cout << "lmat:" << std::endl << lmat << std::endl;

  solver.finalize();
  MPI_Finalize();
  return 0;
}
