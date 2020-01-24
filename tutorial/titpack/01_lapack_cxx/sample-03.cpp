/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// C++ version of TITPACK Ver.2 by H. Nishimori

/************* Sample main program #3 *****************
* 1d Heisenberg antiferromagnet with 16 spins
* Eigenvectors of excited states by lncv1
******************************************************/

#include "titpack.hpp"

int main() {
  std::cout.precision(10);

  // lattice structure
  int n = 16;
  int ibond = n;
  std::vector<int> ipair;
  for (int i = 0; i < ibond; ++i) {
    ipair.emplace_back(i);
    ipair.emplace_back((i + 1) % n);
  }

  // Hamiltonian parameters
  std::vector<double> bondwt(ibond, -1);
  std::vector<double> zrtio(ibond, 1);

  // table of configurations
  std::vector<int> list1;
  std::vector<std::pair<int, int> > list2;
  int idim = sz(n, 0, list1, list2);
  // You may alternatively use szdy or sztn for faster processing
  //   int idim = szdy(n, 0, list1, list2);
  // or
  //   int idim = sztn(n, 0, list1, list2);

  // Eigenvalues
  int nvec = 3;
  int iv = idim / 3 - 1;
  std::vector<double> E, alpha, beta;
  matrix_type coeff;
  matrix_type v;
  int itr = lnc1(n, ipair, bondwt, zrtio, nvec, iv, E, alpha, beta, coeff, v, list1, list2);
  
  std::cout << "[Eigenvalues]\n";
  for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
  std::cout << std::endl;
  std::cout << "[Iteration number]\n\t" << itr << std::endl;

  // Ground-state eigenvector
  matrix_type x;
  lncv1(n, ipair, bondwt, zrtio, nvec, iv, alpha, beta, coeff, x, itr, v, list1, list2);
  // You may alternatively use inv1 / Note: dimension v(4, idim) -
  //   call inv1(n, ipair, bondwt, zrtio, 1, iv, alpha, beta, coeff, x, itr, v, list1, list2);

  std::cout << "[Eigenvector components (selected)]";
  int count = 0;
  for (int i = 12; i < idim; i += idim/20, ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << x(i, nvec - 1);
  }
  std::cout << std::endl;

  // Precision check and correlation functions
  double Hexpec = check1(n, ipair, bondwt, zrtio, x, nvec - 1, v, 0, list1, list2);
}
