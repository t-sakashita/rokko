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

#ifndef TITPACK_MEDIUM_HPP
#define TITPACK_MEDIUM_HPP

#include "common.hpp"
#include "hamiltonian.hpp"

template<typename CRS_MATRIX>
void elm2_mpi(hamiltonian const& hop, CRS_MATRIX& mat) {
  int ibond = hop.num_bonds();
  std::vector<int> cols;
  std::vector<double> values;
  for (int i = mat.start_row(); i <= mat.end_row(); ++i) {
    cols.clear();
    values.clear();
    double diag = 0;

    // diagonal element
    for (int k = 0; k < ibond; ++k) {
      int isite1, isite2;
      boost::tie(isite1, isite2) = hop.site_pair(k);
      int is = (1 << isite1) + (1 << isite2);
      double wght = -hop.bond_weight(k) * hop.z_ratio(k) * 0.5;
      int ibit = hop.config(i) & is;
      diag += (ibit == 0 || ibit == is) ? wght : -wght;
    }
    cols.push_back(i);
    values.push_back(diag);

    // off-diagonal elements
    for (int k = 0; k < ibond; ++k) {
      int isite1, isite2;
      boost::tie(isite1, isite2) = hop.site_pair(k);
      int is = (1 << isite1) + (1 << isite2);
      double wght = -hop.bond_weight(k);
      int ibit = hop.config(i) & is;
      if (ibit != 0 && ibit != is) {
        cols.push_back(hop.config2index(hop.config(i) ^ is));
        values.push_back(wght);
      }
    }
    mat.insert(i, cols, values);
  }
  mat.complete();
}

#endif
