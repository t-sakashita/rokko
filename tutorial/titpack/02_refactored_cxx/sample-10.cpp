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
#include "options.hpp"

int main(int argc, char** argv) {
  std::cout.precision(10);
  options opt(argc, argv, 15);
  if (!opt.valid) std::abort();

  // lattice structure
  int n = opt.N;
  int ibond = n;
  std::vector<int> ipair;
  for (int i = 0; i < ibond; ++i) {
    ipair.emplace_back(i);
    ipair.emplace_back((i + 1) % n);
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
  matrix_type v(hop.dimension(), 3);
  std::vector<double> E, alpha, beta;
  matrix_type coeff;
  matrix_type wk;
  matrix_type x;
  for (int k = 0; k < 3; ++k) {
    int iv = 90 + (hop.dimension() / 3) * k;
    int itr = lnc2(mat, nvec, iv, E, alpha, beta, coeff, wk);

    std::cout << "# " << k << " [Eigenvalues]\n";
    for (int i = 0; i < 4; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl;

    lncv2(mat, nvec, iv, alpha, beta, coeff, x, itr, wk);
    std::cout << "[Eigenvector components (selected)]";
    int count = 0;
    for (int i = 12; i < hop.dimension(); i += hop.dimension()/20, ++count) {
      if (count % 4 == 0) std::cout << std::endl;
      std::cout << '\t' << x(i, 0);
    }
    std::cout << std::endl;
    check2(mat, x, 0, wk, 0);

    for (int i = 0; i < hop.dimension(); ++i) v(i, k) = x(i, 0);
  }

  // Degeneracy check
  std::vector<double> norm;
  int idgn = orthg(v, norm, 3);
  std::cout << " [Degeneracy] : " << idgn << std::endl;
}
