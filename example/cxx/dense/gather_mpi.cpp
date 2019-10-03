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
  std::string solver_name(rokko::parallel_dense_ev::default_solver());
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::grid g;
  if (g.get_myrank() == 0) std::cout << "dimension = " << dim << std::endl;
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());

  rokko::parallel_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::frank_matrix::generate(mat);
  mat.print();

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<matrix_major>> lmat(dim, dim);
  rokko::gather(mat, lmat, 0);
  if (g.get_myrank() == 0)
    std::cout << "lmat:" << std::endl << lmat << std::endl;

  solver.finalize();
  MPI_Finalize();
  return 0;
}
