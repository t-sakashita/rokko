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

#include "common.hpp"
#include <lapacke.h>
#include <iostream>

int sz(int n, double szval, std::vector<int>& list1, std::vector<std::vector<int> >& list2) {
  if (szval < -1e-13 || szval > (n/2 - 1 + 1e-13)) {
    std::cerr << " #(E01)# Variable szval given to sz out of range\n";
    std::abort();
  }

  // initialization
  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;
  int iupspn = n / 2 + (n % 2) + (int)(szval + 0.001);

  list1.clear();
  list2.resize(2);
  list2[0].resize(1 << ihf);
  list2[1].resize(1 << (n - ihf));

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
      list1.push_back(i);
      int ia = i & irght;
      int ib = (i & ilft) / ihfbit;
      if (ib == ibpatn) {
        list2[0][ia] = ja;
        list2[1][ib] = jb;
      } else {
        ibpatn = ib;
        ja = 0;
        jb = icnt;
        list2[0][ia] = ja;
        list2[1][ib] = jb;
      }
      ++ja;
      ++icnt;
    }
  }
  return icnt;
}

void datack(std::vector<int> const& ipair, int n) {
  int ibond = ipair.size() / 2;
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[k * 2];
    int isite2 = ipair[k * 2 + 1];
    if (isite1 < 0 || isite2 < 0 || isite1 >= n || isite2 >= n) {
      std::cerr << " #(E03)# Incorrect data in ipair\n"
                << "         Location :  " << k * 2 << ", " << k * 2 + 1 << std::endl;
      std::abort();
    }
  }
}

boost::tuple<int, int>
bisec(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
      std::vector<double>& E, int ne, double eps, std::vector<int>& iblock,
      std::vector<int>& isplit, double *w) {
  if (E.size() < ne) E.resize(ne);
  if (iblock.size() < ndim) iblock.resize(ndim);
  if (isplit.size() < ndim) isplit.resize(ndim);
  int m, nsplit;
  int info = LAPACKE_dstebz('I', 'B', ndim, 0, 0, 1, ne, eps, &alpha[0], &beta[0],
                            &m, &nsplit, &w[0], &iblock[0], &isplit[0]);
  for (int i = 0; i < ne; ++i) E[i] = w[i];
  return boost::make_tuple(m, nsplit);
}

void vec12(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
           std::vector<double> const& E, int nvec, matrix_type& z,
           std::vector<int>& iblock, std::vector<int>& isplit, double *w) {
  if (z.size1() < nvec || z.size2() != ndim) z.resize(nvec, ndim);
  for (int i = 0; i < nvec; ++i) w[i] = E[i];
  std::vector<int> ifail(nvec);
  int info = LAPACKE_dstein(LAPACK_COL_MAJOR, ndim, &alpha[0], &beta[0], nvec, w, &iblock[0],
                            &isplit[0], &z(0,0), ndim, &ifail[0]);
}

void xcorr(subspace const& ss, std::vector<int> const& npair, const double *x,
           std::vector<double>& sxx) {
  int nbond = npair.size() / 2;
  for (int k=0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
      std::cerr << " #(W01)# Wrong site number given to xcorr\n";
      return;
    }
    double corr = 0;
    int is = (1 << i1) + (1 << i2);
    for (int j = 0; j < ss.dimension(); ++j) {
      int ibit = ss.config(j) & is;
      if (ibit != 0 && ibit != is) corr += x[j] * x[ss.config2index(ss.config(j) ^ is)];
    }
    sxx[k] = corr / 4;
  }
}

void xcorr(subspace const& ss, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& sxx) {
  xcorr(ss, npair, &x[0], sxx);
}

void xcorr(subspace const& ss, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& sxx) {
  xcorr(ss, npair, &x(xindex, 0), sxx);
}

void zcorr(subspace const& ss, std::vector<int> const& npair, const double *x,
           std::vector<double>& szz) {
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
    for (int j = 0; j < ss.dimension(); ++j) {
      int ibit = ss.config(j) & is;
      corr += ((ibit == 0 || ibit == is) ? 1.0 : -1.0) * x[j] * x[j];
    }
    szz[k] = corr / 4;
  }
}

void zcorr(subspace const& ss, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& szz) {
  zcorr(ss, npair, &x(xindex, 0), szz);
}

void zcorr(subspace const& ss, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& szz) {
  zcorr(ss, npair, &x[0], szz);
}

int orthg(matrix_type& ev, std::vector<double>& norm, int numvec) {
  if (numvec <= 1) {
    std::cerr << " #(W03)# Number of vectors is less than 2 in orthg\n";
    return -1;
  }
  if (norm.size() < numvec) norm.resize(numvec);
  int idim = ev.size2();
  for (int i = 0; i < numvec; ++i) {
    double dnorm = 0;
    for (int j = 0; j < idim; ++j) dnorm += ev(i, j) * ev(i, j);
    if (dnorm < 1e-20) {
      std::cerr << " #(W04)# Null vector given to orthg. Location is " << i << std::endl;
      return -1;
    }
    dnorm = 1 / std::sqrt(dnorm);
    for (int j = 0; j < idim; ++j) ev(i, j) *= dnorm;
  }
  int idgn = numvec;
  norm[0] = 1;

  // orthogonalization
  for (int i = 1; i < numvec; ++i) {
    norm[i] = 1;
    for (int j = 0; j < i; ++j) {
      double prjct = 0;
      for (int l = 0; l < idim; ++l) prjct += ev(i, l) * ev(j, l);
      for (int l = 0; l < idim; ++l) ev(i, l) -= prjct * ev(j, l);
    }
    double vnorm = 0;
    for (int l = 0; l < idim; ++l) vnorm += ev(i, l) * ev(i, l);
    if (vnorm > 1e-15) {
      vnorm = 1 / std::sqrt(vnorm);
      for (int l = 0; l < idim; ++l) ev(i, l) *= vnorm;
    } else {
      for (int l = 0; l < idim; ++l) ev(i, l) = 0;
      --idgn;
      norm[i] = 0;
    }
  }

  // check orthogonality
  for (int i = 1; i < numvec; ++i) {
    for (int j = 0; j < i; ++j) {
      double prd = 0;
      for (int l = 0; l < idim; ++l) prd += ev(i, l) * ev(j, l);
      if (std::abs(prd) > 1e-10) {
        std::cerr << " #(W05)# Non-orthogonal vectors at " << i << ' ' << j << std::endl
                  << "         Overlap : " << prd << std::endl;
      }
    }
  }
  return idgn;
}
