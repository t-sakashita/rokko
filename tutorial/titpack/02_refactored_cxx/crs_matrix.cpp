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

#include "crs_matrix.hpp"
#include <tuple>

crs_matrix::crs_matrix(hamiltonian const& hop) {
  int ibond = hop.num_bonds();
  int ic = ibond + 1;

  // initialization
  if (elemnt_.size1() != hop.dimension() || elemnt_.size2() < ic)
    elemnt_.resize(hop.dimension(), ic);
  if (loc_.size1() != hop.dimension() || loc_.size2() < ic) loc_.resize(hop.dimension(), ic);
  for (int i = 0; i < ic; ++i)
    for (int j = 0; j < hop.dimension(); ++j)
      elemnt_(j, i) = 0;

  // diagonal elements
  for (int k = 0; k < ibond; ++k) {
    int isite1, isite2;
    boost::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = -hop.bond_weight(k) * hop.z_ratio(k) * 0.5;
    for (int i = 0; i < hop.dimension(); ++i) {
      int ibit = hop.config(i) & is;
      elemnt_(i, ic - 1) += (ibit == 0 || ibit == is) ? wght : -wght;
      loc_(i, ic - 1) = i;
    }
  }
  
  // off-diagonal elements
  for (int k = 0; k < ibond; ++k) {
    int isite1, isite2;
    boost::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = -hop.bond_weight(k);
    for (int i = 0; i < hop.dimension(); ++i) {
      int ibit = hop.config(i) & is;
      if (ibit == 0 || ibit == is) {
        elemnt_(i, k) = 0;
        loc_(i, k) = i;
      } else {
        elemnt_(i, k) = wght;
        loc_(i, k) = hop.config2index(hop.config(i) ^ is);
      }
    }
  }
}

double crs_matrix::multiply(const double *v1, double *v0) const {
  double prdct = 0;
  for (int k = 0; k < elemnt_.size2(); ++k) {
    for (int j = 0; j < dimension(); ++j) {
      double temp = elemnt_(j, k) * v1[loc_(j, k)];
      v0[j] += temp;
      prdct += v1[j] * temp;
    }
  }
  return prdct;
}
