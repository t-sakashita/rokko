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

#ifndef TITPACK_SMALL_HPP
#define TITPACK_SMALL_HPP

#include "common.hpp"
#include "hamiltonian.hpp"
#include <vector>
#include <boost/tuple/tuple.hpp>

//
// matrix elements
//

template<typename MATRIX>
void elm3(hamiltonian const& hop, MATRIX& elemnt) {
  // initialization
  elemnt.set_zeros();

  // elments
  for (int k = 0; k < hop.num_bonds(); ++k) {
    int isite1, isite2;
    boost::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = hop.bond_weight(k);
    double diag = 0.5 * wght * hop.z_ratio(k);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hop.dimension(); ++i) {
      if (elemnt.is_gindex_myrow(i)) {
        int ibit = hop.config(i) & is;
        if (ibit == 0 || ibit == is) {
          elemnt.update_global(i, i, -diag);
        } else {
          elemnt.update_global(i, i, +diag);
          int newcfg = hop.config2index(hop.config(i) ^ is);
          elemnt.update_global(i, newcfg, -wght);
        }
      }
    }
  }
}

//
// check of the eigenvector and eigenvalue
//
// elemnt    @ nonzero elements
// x         @ eigenvector to be checked
// xindex
// return value: Hexpec <x*H*x>

template<typename MATRIX>
double check3(MATRIX const& elemnt, MATRIX const& x, int xindex) {
  if (elemnt.rows() != elemnt.cols() || elemnt.rows() != x.cols() || x.rows() <= xindex) {
    std::cerr << "check3: Incorrect matrix size\n";
    return 0;
  }
  int idim = elemnt.rows();

  double dnorm = 0;
  for (int j=0; j < idim; ++j) {
    // dnorm += x(xindex, j) * x(xindex, j);
    dnorm += x(j, xindex) * x(j, xindex);
  }
  if (dnorm < 1e-30) {
    std::cerr << " #(W18)# Null vector given to check3\n";
    return 0;
  }
  std::vector<double> v(idim, 0);
      
  for (int j = 0; j < idim; ++j)
    for (int i = 0; i < idim; ++i)
      v[j] += elemnt(i, j) * x(i, xindex);

  double prd = 0;
  // for (int i = 0; i < idim; ++i) prd += v[i] * x(xindex, i);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < idim; ++i) prd += v[i] * x(i, xindex);
  
  std::cout << "---------------------------- Information from check3\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v[i]/x(i, xindex);
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}

#endif
