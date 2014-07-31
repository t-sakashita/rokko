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

/************* Sample main program #10 *****************
* 1d Heisenberg antiferromagnet with 15 spins
* Degeneracy check by various initial vectors
******************************************************/

#include "titpack.hpp"

int main() {
  std::cout.precision(10);

  // lattice structure
  int n = 15;
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
  int idim = sz(n, 0.5, list1, list2);
  // You may alternatively use szdy or sztn for faster processing
  //   int idim = szdy(n, 0, list1, list2);
  // or
  //   int idim = sztn(n, 0, list1, list2);

  // matrix elements
  matrix_type elemnt;
  i_matrix_type loc;
  elm2(n, ipair, bondwt, zrtio, elemnt, loc, list1, list2);

  // Eigenvalues
  int nvec = 1;
  matrix_type v(3, idim);
  for (int k = 0; k < 3; ++k) {
    int iv = 90 + (idim / 3) * k;
    std::vector<double> E, alpha, beta;
    matrix_type coeff;
    matrix_type wk(2, idim);
    int itr = lnc2(elemnt, loc, nvec, iv, E, alpha, beta, coeff, wk);
    
    std::cout << "# " << k << " [Eigenvalues]\n";
    for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl;

    matrix_type x;
    lncv2(elemnt, loc, nvec, iv, alpha, beta, coeff, x, itr, wk);
    std::cout << "[Eigenvector components (selected)]";
    int count = 0;
    for (int i = 12; i < idim; i += idim/20, ++count) {
      if (count % 4 == 0) std::cout << std::endl;
      std::cout << '\t' << x(0, i);
    }
    std::cout << std::endl;
    check2(elemnt, loc, x, 0, wk, 0);

    for (int i = 0; i < idim; ++i) v(k, i) = x(0, i);
  }
  
  // Degeneracy check
  std::vector<double> norm;
  int idgn = orthg(v, norm, 3);
  std::cout << " [Degeneracy] : " << idgn << std::endl;
}
