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

#include "subspace.hpp"
#include <iostream>
#include <cstdlib>

subspace::subspace(int n, double szval) : n_(n), list1_(), list2_() {
  if (szval < -1e-13 || szval > (n/2 - 1 + 1e-13)) {
    std::cerr << " #(E01)# Variable szval given to sz out of range\n";
    std::abort();
  }
  
  // initialization
  int ihf = (n + 1) / 2;
  ihfbit_ = 1 << ihf;
  irght_ = ihfbit_ - 1;
  ilft_ = ((1 << n) - 1) ^ irght_;
  int iupspn = n / 2 + (n % 2) + (int)(szval + 0.001);
  
  list2_.resize(std::max(1 << ihf, 1 << (n - ihf)));

  // main loop
  int icnt = 0;
  int ja = 0;
  int jb = 0;
  int ibpatn = 0;
  for (int i = 0; i < (1 << n); ++i) {
    int isz = 0;
    for (int j = 0; j < n; ++j) {
      isz += (i >> j) & 1;
    }
    if (isz == iupspn) {
      list1_.emplace_back(i);
      int ia = i & irght_;
      int ib = (i & ilft_) / ihfbit_;
      if (ib == ibpatn) {
        list2_[ia].first = ja;
        list2_[ib].second = jb;
      } else {
        ibpatn = ib;
        ja = 0;
        jb = icnt;
        list2_[ia].first = ja;
        list2_[ib].second = jb;
      }
      ++ja;
      ++icnt;
    }
  }
}
