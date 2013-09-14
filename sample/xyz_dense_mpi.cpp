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
#include <fstream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/collective.hpp>

#include <boost/tuple/tuple.hpp>

#include <rokko/utility/xyz_hamiltonian.hpp>
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

  std::ifstream ifs(argv[2]); //str);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    exit(1);
  }

  int L, num_bonds;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > coupling;
  ifs >> L >> num_bonds;
  for (int i=0; i<num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  
  for (int i=0; i<num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(boost::make_tuple(jx, jy, jz));
  }
  
  std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
  for (int i=0; i<num_bonds; ++i) {
    std::cout << lattice[i].first << " " << lattice[i].second << " " << coupling[i].get<0>() << " " << coupling[i].get<1>() << " " << coupling[i].get<2>() << std::endl;
  }
  int dim = 1 << L;
  std::cout << "dim=" << dim << std::endl;
  rokko::solver solver(solver_name);
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  int myrank = g.get_myrank();

  const int root = 0;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
  mat.print();

  rokko::localized_matrix<matrix_major> lmat(dim, dim);
  rokko::gather(mat, lmat, root);
  if (myrank == root)
    std::cout << "lmat:" << std::endl << lmat << std::endl;


  rokko::localized_vector w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g, solver);

  try {
    solver.diagonalize(mat, w, Z);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  // gather of eigenvectors
  rokko::localized_matrix<matrix_major> eigvec_global;
  rokko::localized_matrix<matrix_major> eigvec_sorted(dim, dim);
  rokko::localized_vector eigval_sorted(dim);
  rokko::gather(Z, eigvec_global, root);
  Z.print();
  if (myrank == root) {
    std::cout << "eigvec:" << std::endl << eigvec_global << std::endl;
  }

  std::cout.precision(20);
  /*
  std::cout << "w=" << std::endl;
  for (int i=0; i<dim; ++i) {
    std::cout << w[i] << " ";
  }
  std::cout << std::endl;
  */

  if (myrank == root) {
    rokko::sort_eigenpairs(w, eigvec_global, eigval_sorted, eigvec_sorted);
    std::cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << std::endl;

    /*
    std::cout.precision(3);
    std::cout << "Check the orthogonality of Eigenvectors:" << std::endl
	 << eigvec_sorted * eigvec_sorted.transpose() << std::endl;   // Is it equal to indentity matrix?
    //<< eigvec_global.transpose() * eigvec_global << std::endl;   // Is it equal to indentity matrix?

    std::cout << "residual := A x - lambda x = " << std::endl
         << lmat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << std::endl;
    std::cout << "Are all the following values equal to some eigenvalue = " << std::endl
	 << (lmat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << std::endl;
    //cout << "lmat=" << std::endl << lmat << std::endl;
    */
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
