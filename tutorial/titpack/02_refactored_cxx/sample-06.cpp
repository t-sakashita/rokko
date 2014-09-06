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

/************* Sample main program #6 *****************
* 1d Heisenberg antiferromagnet with 16 spins
* Eigenvalues and an eigenvector / lnc2, lncv2
* Precision check and correlation functions
******************************************************/

#include "titpack.hpp"
#include "options.hpp"

int main(int argc, char** argv) {
  std::cout.precision(10);
  options opt(argc, argv, 16);
  if (!opt.valid) std::abort();

  // lattice structure
  int n = opt.N;
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

  // matrix elements
  crs_matrix mat(hop);
  
  // Eigenvalues
  int nvec = 1;
  int iv = ss.dimension() / 5 - 1;
  std::vector<double> E, alpha, beta;
  matrix_type coeff;
  matrix_type v;
  int itr = lnc2(mat, nvec, iv, E, alpha, beta, coeff, v);
  
  std::cout << "[Eigenvalues]\n";
  for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
  std::cout << std::endl;
  std::cout << "[Iteration number]\n\t" << itr << std::endl;

  // Ground-state eigenvector
  matrix_type x;
  lncv2(mat, 1, iv, alpha, beta, coeff, x, itr, v);

  std::cout << "[Eigenvector components (selected)]";
  int count = 0;
  for (int i = 12; i < ss.dimension(); i += ss.dimension()/20, ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << x(0, i);
  }
  std::cout << std::endl;

  // Precision check and correlation functions
  double Hexpec = check2(mat, x, 0, v, 0);

  std::vector<int> npair;
  npair.push_back(1);
  npair.push_back(2);
  std::vector<double> sxx(1), szz(1);
  xcorr(ss, npair, x, 0, sxx);
  zcorr(ss, npair, x, 0, szz);
  std::cout << "[Nearest neighbor correlation functions]\n\t" 
            << "sxx : " << sxx[0]
            << ", szz : " << szz[0] << std::endl;
}
