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

/************* Sample main program #5 *****************
* 1d Heisenberg antiferromagnet with 16 spins
* Eigenvector of an excited state by inv1
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
  std::vector<std::pair<int, int>> list2;
  int idim = sz(n, 0, list1, list2);
  // You may alternatively use szdy or sztn for faster processing
  //   int idim = szdy(n, 0, list1, list2);
  // or
  //   int idim = sztn(n, 0, list1, list2);

  // Eigenvalues
  int nvec = 1;
  matrix_type v(idim, 2);
  std::vector<double> E, alpha, beta;
  matrix_type coeff;
  matrix_type wk;
  matrix_type x;
  for (int k = 0; k < 2; ++k) {
    int iv = 20 + (idim / 2) * k;
    int itr = lnc1(n, ipair, bondwt, zrtio, nvec, iv, E, alpha, beta, coeff, wk, list1, list2);
    
    std::cout << "# " << k << " [Eigenvalues]\n";
    for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl;

    lncv1(n, ipair, bondwt, zrtio, nvec, iv, alpha, beta, coeff, x, itr, wk, list1, list2);
    for (int i = 0; i < idim; ++i) v(i, k) = x(i, 0);
  }
  
  // Degeneracy check
  std::vector<double> norm;
  int idgn = orthg(v, norm, 2);
  std::cout << " [Degeneracy] : " << idgn << "   Norm : " << norm[0] << '\t' << norm[1] << std::endl;
}
