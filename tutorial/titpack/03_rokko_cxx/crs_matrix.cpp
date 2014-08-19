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
#include <boost/tuple/tuple.hpp>

crs_matrix::crs_matrix(hamiltonian const& hop) {
  int ibond = hop.num_bonds();
  int ic = ibond + 1;

  // initialization
  if (elemnt_.rows() < ic || elemnt_.cols() != hop.dimension())
    elemnt_.resize(ic, hop.dimension());
  if (loc_.size1() < ic || loc_.size2() != hop.dimension()) loc_.resize(ic, hop.dimension());
  for (int i = 0; i < ic; ++i)
    for (int j = 0; j < hop.dimension(); ++j)
      elemnt_(i, j) = 0;

  // diagonal elements
  for (int k = 0; k < ibond; ++k) {
    int isite1, isite2;
    boost::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = -hop.bond_weight(k) * hop.z_ratio(k) * 0.5;
    for (int i = 0; i < hop.dimension(); ++i) {
      int ibit = hop.config(i) & is;
      elemnt_(ic - 1, i) += (ibit == 0 || ibit == is) ? wght : -wght;
      loc_(ic - 1, i) = i;
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
        elemnt_(k, i) = 0;
        loc_(k, i) = i;
      } else {
        elemnt_(k, i) = wght;
        loc_(k, i) = hop.config2index(hop.config(i) ^ is);
      }
    }
  }
}

double crs_matrix::multiply(const double *v1, double *v0) const {
  double prdct = 0;
  for (int k = 0; k < elemnt_.rows(); ++k) {
    for (int j = 0; j < dimension(); ++j) {
      double temp = elemnt_(k, j) * v1[loc_(k, j)];
      v0[j] += temp;
      prdct += v1[j] * temp;
    }
  }
  return prdct;
}
