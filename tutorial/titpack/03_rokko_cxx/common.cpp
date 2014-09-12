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

#include "common.hpp"
#include <lapacke.h>
#include <iostream>

void xcorr(subspace const& ss, std::vector<int> const& npair, const double *x,
           std::vector<double>& sxx) {
  int nbond = npair.size() / 2;
  for (int k = 0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
      std::cerr << " #(W01)# Wrong site number given to xcorr\n";
      return;
    }
    double corr = 0;
    int is = (1 << i1) + (1 << i2);
    for (int j = 0; j < ss.dimension(); ++j) {
      int ibit = ss.config(j) & is;
      if (ibit != 0 && ibit != is) corr += x[j] * x[ss.config2index(ss.config(j) ^ is)];
    }
    sxx[k] = corr / 4;
  }
}

void xcorr(subspace const& ss, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& sxx) {
  xcorr(ss, npair, &x[0], sxx);
}
