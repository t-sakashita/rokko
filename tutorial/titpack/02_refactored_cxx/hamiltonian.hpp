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

#ifndef TITPACK_HAMILTONIAN_HPP
#define TITPACK_HAMILTONIAN_HPP

#include "subspace.hpp"

//
// hamiltonian
//

class hamiltonian {
public:
  hamiltonian(subspace const& ss, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio);
  int num_sites() const { return ss_.num_sites(); }
  int num_bonds() const { return bondwt_.size(); }
  int dimension() const { return ss_.dimension(); }
  int config(int i) const { return ss_.config(i); }
  int config2index(int c) const { return ss_.config2index(c); }
  std::pair<int, int> site_pair(int k) const {
    return std::make_pair(ipair_[2 * k], ipair_[2 * k + 1]);
  }
  double bond_weight (int k) const { return bondwt_[k]; }
  double z_ratio (int k) const { return zrtio_[k]; }
  double multiply(const double *v1, double *v0) const;
private:
  subspace const& ss_;
  std::vector<int> ipair_;
  std::vector<double> bondwt_, zrtio_;
};

#endif
