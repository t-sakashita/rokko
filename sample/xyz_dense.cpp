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

#include <iostream>
#include <fstream>

#include <rokko/serial_solver.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/collective.hpp>

#include <boost/tuple/tuple.hpp>

#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

#include <boost/lexical_cast.hpp>

int main(int argc, char *argv[]) {
  //typedef rokko::matrix_row_major matrix_major;
  typedef rokko::matrix_col_major matrix_major;

  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name matrix_size" << std::endl;
    exit(34);
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
  rokko::serial_solver solver(solver_name);
  solver.initialize(argc, argv);

  rokko::localized_matrix<rokko::matrix_col_major> mat(dim, dim);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
  std::cout << mat << std::endl;

  rokko::localized_vector w(dim);
  rokko::localized_matrix<rokko::matrix_col_major> eigvec(dim, dim);

  try {
    solver.diagonalize(mat, w, eigvec);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }

  std::cout << "eigvec:" << eigvec << std::endl;
  
  std::cout.precision(20);
  rokko::localized_vector eigval_sorted(dim);
  rokko::localized_matrix<rokko::matrix_col_major> eigvec_sorted(dim, dim);
  rokko::sort_eigenpairs(w, eigvec, eigval_sorted, eigvec_sorted);
  std::cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << std::endl;
  
  std::cout.precision(3);
  std::cout << "Check the orthogonality of Eigenvectors:" << std::endl
            << eigvec_sorted * eigvec_sorted.transpose() << std::endl;   // Is it equal to indentity matrix?
  //<< eigvecs.transpose() * eigvecs << std::endl;   // Is it equal to indentity matrix?
  
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
  std::cout << "residual := A x - lambda x = " << std::endl
         << mat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << std::endl;
  std::cout << "Are all the following values equal to some eigenvalue = " << std::endl
            << (mat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << std::endl;
  //cout << "lmat=" << std::endl << lmat << std::endl;

  solver.finalize();
  return 0;
}
