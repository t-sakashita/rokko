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

#include <iostream>
#include <fstream>
#include <boost/tuple/tuple.hpp>

#include <rokko/serial_solver.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

int main(int argc, char *argv[]) {
try {
  typedef rokko::matrix_col_major matrix_major;

  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name lattice_file" << std::endl;
    exit(34);
  }
  std::string solver_name(argv[1]);
  std::cout << "solver = " << solver_name << std::endl;
  std::ifstream ifs(argv[2]);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    exit(1);
  }
  std::cout << "lattice file = " << argv[2] << std::endl;

  rokko::serial_solver solver(solver_name);
  solver.initialize(argc, argv);

  int N, num_bonds;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > coupling;
  ifs >> N >> num_bonds;
  for (int i = 0; i < num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  
  for (int i = 0; i < num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(boost::make_tuple(jx, jy, jz));
  }

  std::cout << "number of sites = " << N << " number of bonds = " << num_bonds << std::endl;
  int dim = 1 << N;
  std::cout << "matrix dimension = " << dim << std::endl;
  rokko::localized_matrix<matrix_major> mat(dim, dim);
  rokko::xyz_hamiltonian::generate(N, lattice, coupling, mat);

  rokko::localized_matrix<matrix_major> mat_in(mat);
  rokko::localized_vector eigval(dim);
  rokko::localized_matrix<matrix_major> eigvec(dim, dim);

  solver.diagonalize(mat_in, eigval, eigvec);

  rokko::localized_vector eigval_sorted(dim);
  rokko::localized_matrix<matrix_major> eigvec_sorted(dim, dim);
  rokko::sort_eigenpairs(eigval, eigvec, eigval_sorted, eigvec_sorted);

  std::cout.precision(10);
  std::cout << "10 Lowest Eigenvalues (value, degeneracy, residue)" << std::endl;
  int index = dim - 1;
  int deg = 0;
  double val = eigval_sorted(index);
  double res = 0;
  for (int i = 0; i < 10; ++i) {
    for (; index >= 0; --index) {
      bool change = (std::abs((eigval_sorted(index) - val) / val) > 1e-10);
      if (!change) {
        ++deg;
        res = std::max(res, std::abs(eigvec_sorted.transpose().row(index) * mat *
                                     eigvec_sorted.col(index) - eigval_sorted(index)));
      }
      if (change || index == 0) {
        std::cout << i + 1 << ":\t" << val << '\t' << deg << '\t' << res << std::endl;
        val = eigval_sorted(index);
        deg = 1;
        res = std::abs(eigvec_sorted.transpose().row(index) * mat *
                       eigvec_sorted.col(index) - eigval_sorted(index));
        break;
      }
    }
    if (index == 0) break;
  }

  solver.finalize();
  return 0;
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
}
