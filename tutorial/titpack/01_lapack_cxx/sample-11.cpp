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

/************ Sample main program #11 *****************
* 1d Heisenberg antiferromagnet with 8 spins
* Eigenvalues and an eigenvector by diag
******************************************************/

#include "titpack.hpp"

int main() {
  std::cout.precision(10);

  // lattice structure
  int n = 8;
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
  //   int = szdy(n, 0, list1, list2);
  // or
  //   int idim = sztn(n, 0, list1, list2);

  // Hamiltonian matrix
  matrix_type elemnt;
  elm3(n, ipair, bondwt, zrtio, elemnt, list1, list2);
  
  std::vector<double> E;
  matrix_type v;
  int nvec = 1;
  diag(elemnt, E, v, nvec);

  int ne = 4;
  std::cout << "[Eigenvalues]\n";
  for (int i = 0; i < ne; ++i) std::cout << '\t' << E[i];
  std::cout << std::endl;

  // // Do not forget to call elm3 again before calling check3
  elm3(n, ipair, bondwt, zrtio, elemnt, list1, list2);
  check3(elemnt, v, 0);
  
  std::vector<int> npair;
  npair.push_back(1);
  npair.push_back(2);
  std::vector<double> sxx(1);
  xcorr(n, npair, v, 0, sxx, list1, list2);
  std::cout << "sxx: " << sxx[0] << std::endl;
  std::vector<double> szz(1);
  zcorr(n, npair, v, 0, szz, list1);
  std::cout << "szz: " << szz[0] << std::endl;
}
