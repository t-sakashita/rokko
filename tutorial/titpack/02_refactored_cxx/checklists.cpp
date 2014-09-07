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

#include "common.hpp"

void test(int n, double total_sz) {
  std::vector<int> list1;
  std::vector<std::vector<int> > list2;
  int idim = sz(n, total_sz, list1, list2);
  std::cout << n << ' ' << total_sz << ' ' << idim << std::endl;
  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;
  for (int i = 0; i < list1.size(); ++i) {
    int ia = list1[i] & irght;
    int ib = (list1[i] & ilft) / ihfbit;
    if (i != list2[0][ia] + list2[1][ib]) {
      std::cerr << "Incorrect table list1 and list2\n";
      std::abort();
    }
  }
}

int main() {
  test(16, 0);
  test(8, 0);
  test(8, 1);
}
