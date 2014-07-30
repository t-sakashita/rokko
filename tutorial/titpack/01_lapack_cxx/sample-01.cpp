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

/************* Sample main program #1 *****************
* 1d Heisenberg antiferromagnet with 16 spins
* Eigenvalues and an eigenvector / lnc1, lncv1
* Precision check and correlation functions
******************************************************/

#include "titpack.hpp"

int main() {
  std::cout.precision(10);

  // lattice structure
  int n = 16;
  int ibond = n;
  std::vector<int> ipair;
  for (int i = 0; i < ibond; ++i) {
    ipair.push_back(i);
    ipair.push_back((i + 1) % n);
  }

  // Hamiltonian parameters
  std::vector<double> bondwt(ibond, -1);
  std::vector<double> zrtio(ibond, 1);

  // table of configurations
  std::vector<int> list1;
  std::vector<std::vector<int> > list2;
  int idim = sz(n, 0, list1, list2);
  // You may alternatively use szdy or sztn for faster processing
  //   int idim = szdy(n, 0, list1, list2);
  // or
  //   int idim = sztn(n, 0, list1, list2);

  // Eigenvalues
  int nvec = 1;
  int iv = idim / 3 - 1;
  std::vector<double> E, alpha, beta;
  matrix_type coeff;
  matrix_type v(2, idim);
  int itr = lnc1(n, ipair, bondwt, zrtio, nvec, iv, E, alpha, beta, coeff, v, list1, list2);
  
  std::cout << "[Eigenvalues]\n";
  for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
  std::cout << std::endl;
  std::cout << "[Iteration number]\n\t" << itr << std::endl;

  // Ground-state eigenvector
  matrix_type x;
  lncv1(n, ipair, bondwt, zrtio, 1, iv, alpha, beta, coeff, x, itr, v, list1, list2);
  // You may alternatively use inv1 / Note: dimension v(4, idim) -
  //   call inv1(n, ipair, bondwt, zrtio, 1, iv, alpha, beta, coeff, x, itr, v, list1, list2);

  std::cout << "[Eigenvector components (selected)]";
  int count = 0;
  for (int i = 12; i < idim; i += idim/20, ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << x(0, i);
  }
  std::cout << std::endl;

  // Precision check and correlation functions
  double Hexpec = check1(n, ipair, bondwt, zrtio, x, 0, v, 0, list1, list2);

  std::vector<int> npair;
  npair.push_back(1);
  npair.push_back(2);
  std::vector<double> sxx(1), szz(1);
  xcorr(n, npair, x, 0, sxx, list1, list2);
  zcorr(n, npair, x, 0, szz, list1);
  std::cout << "[Nearest neighbor correlation functions]\n\t" 
            << "sxx : " << sxx[0]
            << ", szz : " << szz[0] << std::endl;
}
