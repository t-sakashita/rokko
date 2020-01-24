/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
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
#include "options.hpp"
#include <chrono>

int main(int argc, char** argv) {
  std::cout.precision(10);
  options opt(argc, argv, 14);
  if (!opt.valid) std::abort();
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

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

  // Hamiltonian matrix
  matrix_type elemnt;
  elm3(hop, elemnt);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  
  std::vector<double> E;
  matrix_type v;
  int nvec = 1;
  diag(elemnt, E, v, nvec);
  std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();

  int ne = 4;
  std::cout << "[Eigenvalues]\n";
  for (int i = 0; i < ne; ++i) std::cout << '\t' << E[i];
  std::cout << std::endl;

  // // Do not forget to call elm3 again before calling check3
  elm3(hop, elemnt);
  check3(elemnt, v, 0);
  std::chrono::system_clock::time_point t4 = std::chrono::system_clock::now();
  
  std::vector<int> npair;
  npair.emplace_back(1);
  npair.emplace_back(2);
  std::vector<double> sxx(1);
  xcorr(ss, npair, v, 0, sxx);
  std::cout << "sxx: " << sxx[0] << std::endl;
  std::vector<double> szz(1);
  zcorr(ss, npair, v, 0, szz);
  std::cout << "szz: " << szz[0] << std::endl;
  std::chrono::system_clock::time_point t5 = std::chrono::system_clock::now();
  std::cerr << "initialize      " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << " sec\n"
            << "diagonalization " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count() << " sec\n"
            << "check           " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count() << " sec\n"
            << "correlation     " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t5-t4).count() << " sec\n";
}
