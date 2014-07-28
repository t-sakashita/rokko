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

#include "small.hpp"
#include "common.hpp"
#include <lapacke.h>
#include <iostream>

void elm3(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, matrix_type& elemnt,
          std::vector<int> const& list1, std::vector<std::vector<int> > const& list2) {
  int ihfbit = 1 << ((n+1)/2);
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;
  int idim = list1.size();

  int ibond = ipair.size() / 2;
  if (ipair.size() != 2 * ibond || bondwt.size() != ibond || zrtio.size() != ibond) {
    std::cerr << "Incorrect size of ipair, bondwt, or zrtio\n";
    return;
  }
  datack(ipair, n);
  
  // initialization
  for (int i = 0; i < idim; ++i)
    for (int j = 0; j < idim; ++j)
      elemnt(i, j) = 0;

  // elments
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[2 * k];
    int isite2 = ipair[2 * k + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = bondwt[k];
    double diag = 0.5 * wght * zrtio[k];
    for (int i = 0; i < idim; ++i) {
      int ibit = list1[i] & is;
      if (ibit == 0 || ibit == is) {
        elemnt(i,i) -= diag;
      } else {
        elemnt(i,i) += diag;
        int iexchg = list1[i] ^ is;
        int ia = iexchg & irght;
        int ib = (iexchg & ilft) / ihfbit;
        int newcfg = list2[0][ia] + list2[1][ib];
        elemnt(i, newcfg) = -wght;
      }
    }
  }
}

void diag(matrix_type& elemnt, std::vector<double>& E, matrix_type& v, int nvec) {
  if (elemnt.size1() != elemnt.size2()) {
    std::cerr << "diag: Incorrect matrix size\n";
    return;
  }
  int idim = elemnt.size1();
  E.resize(idim);
  int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', idim, &elemnt(0,0), idim, &E[0]);

  v.resize(nvec, idim);
  for (int j = 0; j < nvec; ++j)
    for (int i = 0; i < idim; ++i)
      v(j, i) = elemnt(j, i);
}

double check3(matrix_type const& elemnt, matrix_type const& x, int xindex) {
  if (elemnt.size1() != elemnt.size2() || elemnt.size1() != x.size2() || x.size1() <= xindex) {
    std::cerr << "check3: Incorrect matrix size\n";
    return 0;
  }
  int idim = elemnt.size1();

  double dnorm = 0;
  for (int j=0; j < idim; ++j) {
    dnorm += x(xindex, j) * x(xindex, j);
  }
  if (dnorm < 1e-30) {
    std::cerr << " #(W18)# Null vector given to check3\n";
    return 0;
  }
  std::vector<double> v(idim, 0);
      
  for (int j = 0; j < idim; ++j)
    for (int i = 0; i < idim; ++i)
      v[j] += elemnt(j, i) * x(xindex, i);

  double prd = 0;
  for (int i = 0; i < idim; ++i) prd += v[i] * x(xindex, i);
  
  std::cout << "---------------------------- Information from check3\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13),idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v[i]/x(xindex, i);
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}
