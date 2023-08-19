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

auto magnetization(int L, std::vector<std::pair<int, int>>& lattice, int power, const double* v) {
  const int N = 1 << L;
  const int m_power = 1 << power;
  const double coeff = 1. / m_power;
  double sum = 0;

  for (int k=0; k<N; ++k) {
    double m_z = 0;
    for (int i=0; i<L; ++i) {
      const int mask = 1 << i;
      if ((k & mask) == mask) {
        m_z += 0.5;
      }
      else {
        m_z -=0.5;
      }
    }      
    sum += std::pow(m_z, static_cast<double>(power)) * v[k] * v[k];
  }

  return sum;
}

auto magnetization(int L, std::vector<std::pair<int, int>>& lattice, int power, const std::vector<double>& v) {
  return magnetization(L, lattice, power, v.data());
}

//} // namespace heisenberg_hamiltonian

} // namespace rokko
