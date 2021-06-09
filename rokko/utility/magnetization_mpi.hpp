/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <vector>
#include <rokko/eigen3.hpp>

namespace rokko {

//namespace heisenberg_hamiltonian {

double magnetization(const MPI_Comm& comm, int L, std::vector<std::pair<int, int>>& lattice, int power, const double* v) {
  int myrank, nproc;
  MPI_Status status;
  int ierr;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  int n = nproc;
  int r = myrank;
  int p = -1;
  int p_count = 0;
  do {
    n /= 2;
    ++p;
    r = r >> 1;
    p_count += r & 1; 
  } while (n > 0);

  std::cout << "myrank=" << myrank << "p_count=" << p_count << std::endl;
 
 int N = 1 << L;
  int  m_power = 1 << power;
  double coeff = 1. / m_power;
  double sum = 0;

  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int mask = 1 << i;
    for (int k=0; k<N; ++k) {
      if ((k & mask) == mask) {
        sum += coeff * v[k] * v[k];
      }
      else {
        sum += - coeff * v[k] * v[k];
      }
    }
  }

  int i = lattice[lattice.size()-1].second;
  int mask = 1 << i;
  for (int k=0; k<N; ++k) {
    if ((k & mask) == mask) {
      sum += coeff * v[k] * v[k];
    }
    else {
      sum += - coeff * v[k] * v[k];
    }
  }

  return sum;
}

double magnetization(int L, std::vector<std::pair<int, int>>& lattice, int power, const std::vector<double>& v) {
  return magnetization(L, lattice, power, v.data());
}

//} // namespace heisenberg_hamiltonian

} // namespace rokko
