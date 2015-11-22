/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <fstream>
#include <rokko/rokko.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#define BOOST_TEST_MODULE test_solver
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_solver) {
  std::vector<std::string> names;
  int argc = boost::unit_test::framework::master_test_suite().argc;
  if (argc == 1) {
    names = rokko::serial_dense_ev::solvers();
  } else {
    for (int num=1; num < argc; ++num) {
      names.push_back(boost::unit_test::framework::master_test_suite().argv[num]);
    }
  }

  BOOST_FOREACH(std::string name, names) {
    std::vector<std::pair<int, int> > lattice;
    std::vector<boost::tuple<double, double, double> > coupling;
    int L = 5;
    int num_bonds = L-1;
    int dim = 1 << L;

    for (int i=0; i<num_bonds; ++i) {
      lattice.push_back(std::make_pair(i, i+1));
      coupling.push_back(boost::make_tuple(1, 1, 1));
    }

    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (int i=0; i<num_bonds; ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << " " << coupling[i].get<0>() << " " << coupling[i].get<1>() << " " << coupling[i].get<2>() << std::endl;
    }

    std::cout << "solver=" << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    rokko::localized_matrix<double, rokko::matrix_col_major> mat(dim, dim);
    rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
    rokko::localized_vector<double> w(dim);
    rokko::localized_matrix<double, rokko::matrix_col_major> Z(dim, dim);
    std::cout << "mat=" << mat << std::endl;
    solver.diagonalize(mat, w, Z);
    
    double sum = 0;
    for(int i=0; i<dim; ++i) {
      sum += w[i];
    }

    std::cout << "w=" << w.transpose() << std::endl;
    //BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);
    
    solver.finalize();
  }
}
