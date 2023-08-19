/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
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

using matrix_major = rokko::matrix_col_major;
// using matrix_major = rokko::matrix_row_major;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  unsigned int dim = 10;
  std::string solver_name(rokko::parallel_dense_ev::default_solver());
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) dim = std::stoi(argv[2]);

  rokko::grid g;
  if (g.get_myrank() == 0) std::cout << "dimension = " << dim << std::endl;
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());

  rokko::parallel_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  const rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::frank_matrix::generate(mat);
  mat.print();

  Eigen::MatrixXd lmat(dim, dim);
  for (int proc = 0; proc < g.get_nprocs(); ++proc) {
    rokko::gather(mat, lmat, proc);
    if (g.get_myrank() == proc) {
      std::cout << "root = " << proc << std::endl;
      std::cout << "lmat:" << std::endl << lmat << std::endl;
    }
    std::cout << std::flush;
    MPI_Barrier(g.get_comm());
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
