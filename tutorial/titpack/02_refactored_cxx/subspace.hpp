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

#ifndef TITPACK_SUBSPACE_HPP
#define TITPACK_SUBSPACE_HPP

#include <vector>

//
// subspace : configurations with the specified sz
//

class subspace {
public:
  subspace(int n, double szval);
  int num_sites() const { return n_; }
  int dimension() const { return list1_.size(); }
  int config(int i) const { return list1_[i]; }
  int config2index(int c) const {
    int ia = c & irght_;
    int ib = (c & ilft_) / ihfbit_;
    return list2_[ia].first + list2_[ib].second;
  }
private:
  int n_, ihfbit_, irght_, ilft_;
  std::vector<int> list1_;
  std::vector<std::pair<int, int> > list2_;
};

#endif
