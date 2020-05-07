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

#ifndef TITPACK_COMMON_HPP
#define TITPACK_COMMON_HPP

#include "subspace.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>

//
// xx correlation function
//
// n           @ lattice size
// npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
// x           @ eigenvetor
// sxx         # xx correlation function
// list1,list2 @ spin configurations generated in 'sz'

void xcorr(subspace const& ss, std::vector<int> const& npair, const double *x,
           std::vector<double>& sxx);

void xcorr(subspace const& ss, std::vector<int> const& npair,
           std::vector<double> const& x, std::vector<double>& sxx);

template<typename MATRIX>
void xcorr(subspace const& ss, std::vector<int> const& npair,
           MATRIX const& x, int xindex, std::vector<double>& sxx) {
  xcorr(ss, npair, &x(0, xindex), sxx);
}

//
// ************* zz correlation function **************
//

template<typename MATRIX>
void zcorr(subspace const& ss, std::vector<int> const& npair,
           MATRIX const& x, int xindex, std::vector<double>& szz) {
  int nbond = npair.size() / 2;
  for (int k = 0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
      std::cerr << " #(W02)# Wrong site number given to zcorr\n";
      return;
    }
    double corr = 0;
    int is = (1 << i1) + (1 << i2);
    for (int i = 0; i < ss.dimension(); ++i) {
      int ibit = ss.config(i) & is;
      corr += ((ibit == 0 || ibit == is) ? 1.0 : -1.0) * x(i, xindex) * x(i, xindex);
    }
    szz[k] = corr / 4;
  }
}

template<typename MATRIX>
void zcorr_mpi(subspace const& ss, std::vector<int> const& npair,
               MATRIX const& x, int xindex, std::vector<double>& szz) {
  int nbond = npair.size() / 2;
  std::vector<double> szz_local(nbond);
  if (x.has_global_col_index(xindex)) {
    for (int k = 0; k < nbond; ++k) {
      int i1 = npair[k * 2];
      int i2 = npair[k * 2 + 1];
      if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
        std::cerr << " #(W02)# Wrong site number given to zcorr\n";
        return;
      }
      double corr = 0;
      int is = (1 << i1) + (1 << i2);
      for (int local_i = 0; local_i < x.get_m_local(); ++local_i) {
        int i = x.translate_l2g_row(local_i);
        int ibit = ss.config(i) & is;
        double val = x.get_global(i, xindex);
        corr += ((ibit == 0 || ibit == is) ? 1.0 : -1.0) * val * val;
      }
      szz_local[k] = corr / 4;
    }
  }
  MPI_Reduce(szz_local.data(), szz.data(), nbond, MPI_DOUBLE, MPI_SUM, 0, x.get_grid().get_comm());
}

//
// Orthogonalization of the eigenvectors
//
// return value #  degree of degenearcy
// ev      @# vectors to be orthogonalized / orthogonalized vectors
// norm(j) #  norm of the j-th vector returned
// numvec  @  number of vectors to be checked

template<typename MATRIX>
int orthg(MATRIX& ev, std::vector<double>& norm, int numvec) {
  if (numvec <= 1) {
    std::cerr << " #(W03)# Number of vectors is less than 2 in orthg\n";
    return -1;
  }
  if (norm.size() < numvec) norm.resize(numvec);
  int idim = ev.rows();
  for (int i = 0; i < numvec; ++i) {
    double dnorm = 0;
    for (int j = 0; j < idim; ++j) dnorm += ev(j, i) * ev(j, i);
    if (dnorm < 1e-20) {
      std::cerr << " #(W04)# Null vector given to orthg. Location is " << i << std::endl;
      return -1;
    }
    dnorm = 1 / std::sqrt(dnorm);
    for (int j = 0; j < idim; ++j) ev(j, i) *= dnorm;
  }
  int idgn = numvec;
  norm[0] = 1;

  // orthogonalization
  for (int i = 1; i < numvec; ++i) {
    norm[i] = 1;
    for (int j = 0; j < i; ++j) {
      double prjct = 0;
      for (int l = 0; l < idim; ++l) prjct += ev(l, i) * ev(l, j);
      for (int l = 0; l < idim; ++l) ev(l, i) -= prjct * ev(l, j);
    }
    double vnorm = 0;
    for (int l = 0; l < idim; ++l) vnorm += ev(l, i) * ev(l, i);
    if (vnorm > 1e-15) {
      vnorm = 1 / std::sqrt(vnorm);
      for (int l = 0; l < idim; ++l) ev(l, i) *= vnorm;
    } else {
      for (int l = 0; l < idim; ++l) ev(l, i) = 0;
      --idgn;
      norm[i] = 0;
    }
  }

  // check orthogonality
  for (int i = 1; i < numvec; ++i) {
    for (int j = 0; j < i; ++j) {
      double prd = 0;
      for (int l = 0; l < idim; ++l) prd += ev(l, i) * ev(l, j);
      if (std::abs(prd) > 1e-10) {
        std::cerr << " #(W05)# Non-orthogonal vectors at " << i << ' ' << j << std::endl
                  << "         Overlap : " << prd << std::endl;
      }
    }
  }
  return idgn;
}

#endif
