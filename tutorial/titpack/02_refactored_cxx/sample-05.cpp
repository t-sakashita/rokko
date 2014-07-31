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
    ipair.push_back(i);
    ipair.push_back((i + 1) % n);
  }

  // Hamiltonian parameters
  std::vector<double> bondwt(ibond, -1);
  std::vector<double> zrtio(ibond, 1);

  // table of configurations and Hamiltonian operator
  subspace ss(n, 0);
  hamiltonian hop(ss, ipair, bondwt, zrtio);

  // Eigenvalues
  int nvec = 1;
  matrix_type v(2, ss.dimension());
  for (int k = 0; k < 2; ++k) {
    int iv = 20 + (ss.dimension() / 2) * k;
    std::vector<double> E, alpha, beta;
    matrix_type coeff;
    matrix_type wk;
    int itr = lnc1(hop, nvec, iv, E, alpha, beta, coeff, wk);
    
    std::cout << "# " << k << " [Eigenvalues]\n";
    for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl;

    matrix_type x;
    lncv1(hop, nvec, iv, alpha, beta, coeff, x, itr, wk);
    for (int i = 0; i < ss.dimension(); ++i) v(k, i) = x(0, i);
  }
  
  // Degeneracy check
  std::vector<double> norm;
  int idgn = orthg(v, norm, 2);
  std::cout << " [Degeneracy] : " << idgn << "   Norm : " << norm[0] << '\t' << norm[1] << std::endl;
}
