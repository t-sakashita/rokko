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

#include "hamiltonian.hpp"

hamiltonian::hamiltonian(subspace const& ss, std::vector<int> const& ipair,
                         std::vector<double> const& bondwt, std::vector<double> const& zrtio) :
  ss_(ss), ipair_(ipair), bondwt_(bondwt), zrtio_(zrtio) {}

void hamiltonian::multiply(const double *v1, double *v0) const {
  for (int k = 0; k < num_bonds(); ++k) {
    int isite1 = ipair_[k * 2];
    int isite2 = ipair_[k * 2 + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = 0.5 * bondwt_[k] * zrtio_[k];
    for (int j = 0; j < dimension(); ++j) {
      int ibit = config(j) & is;
      double factor, offdg;
      if (ibit == 0 || ibit == is) {
        v0[j] -= wght * v1[j];
        factor = 1;
        offdg = 0;
      } else {
        v0[j] += wght * v1[j];
        double temp = v1[config2index(config(j) ^ is)] * bondwt_[k];
        v0[j] -= temp;
        factor = -1;
        offdg = -temp * v1[j];
      }
    }
  }
  return;
}
