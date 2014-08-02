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
#include <lapacke.h>
#include <iostream>
#include <boost/tuple/tuple.hpp>

void elm3(hamiltonian const& hop, matrix_type& elemnt) {
  // initialization
  elemnt.resize(hop.dimension(), hop.dimension());
  for (int i = 0; i < hop.dimension(); ++i)
    for (int j = 0; j < hop.dimension(); ++j)
      elemnt(i, j) = 0;

  // elments
  for (int k = 0; k < hop.num_bonds(); ++k) {
    int isite1, isite2;
    boost::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = hop.bond_weight(k);
    double diag = 0.5 * wght * hop.z_ratio(k);
    for (int i = 0; i < hop.dimension(); ++i) {
      int ibit = hop.config(i) & is;
      if (ibit == 0 || ibit == is) {
        elemnt(i,i) -= diag;
      } else {
        elemnt(i,i) += diag;
        int newcfg = hop.config2index(hop.config(i) ^ is);
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
  if (E.size() < idim) E.resize(idim);
  int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', idim, &elemnt(0,0), idim, &E[0]);

  if (v.size1() < nvec || v.size2() != idim) v.resize(nvec, idim);
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
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v[i]/x(xindex, i);
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}
