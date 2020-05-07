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

int sz(int n, double szval, std::vector<int>& list1, std::vector<std::pair<int, int>>& list2) {
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
  list2.resize(1 << std::max(ihf, n - ihf));

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
      list1.emplace_back(i);
      int ia = i & irght;
      int ib = (i & ilft) / ihfbit;
      if (ib == ibpatn) {
        list2[ia].first = ja;
        list2[ib].second = jb;
      } else {
        ibpatn = ib;
        ja = 0;
        jb = icnt;
        list2[ia].first = ja;
        list2[ib].second = jb;
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

std::tuple<int, int>
bisec(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
      std::vector<double>& E, int ne, double eps, std::vector<int>& iblock,
      std::vector<int>& isplit, double *w) {
  if (E.size() < ne) E.resize(ne);
  if (iblock.size() < ndim) iblock.resize(ndim);
  if (isplit.size() < ndim) isplit.resize(ndim);
  int m, nsplit;
  int info = LAPACKE_dstebz('I', 'B', ndim, 0, 0, 1, ne, eps, alpha.data(), beta.data(),
                            &m, &nsplit, w.data(), iblock.data(), isplit.data());
  for (int i = 0; i < ne; ++i) E[i] = w[i];
  return std::make_tuple(m, nsplit);
}

void vec12(std::vector<double> const& alpha, std::vector<double> const& beta, int ndim,
           std::vector<double> const& E, int nvec, matrix_type& z,
           std::vector<int>& iblock, std::vector<int>& isplit, double *w) {
  if (z.size1() != ndim || z.size2() < nvec) z.resize(ndim, nvec);
  for (int i = 0; i < nvec; ++i) w[i] = E[i];
  std::vector<int> ifail(nvec);
  int info = LAPACKE_dstein(LAPACK_COL_MAJOR, ndim, alpha.data(), beta.data(), nvec, w, iblock.data(),
                            isplit.data(), &z(0,0), ndim, ifail.data());
}

void xcorr(int n, std::vector<int> const& npair, const double *x,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2) {
  int idim = list1.size();
  int nbond = npair.size() / 2;
  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = (1 << ihf) - 1;
  int ilft = ((1 << n) - 1) ^ irght;

  for (int k=0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    if (i1 < 0 || i1 >= n || i2 < 0 || i2 >= n || i1 == i2) {
      std::cerr << " #(W01)# Wrong site number given to xcorr\n";
      return;
    }
    double corr = 0;
    int is = (1 << i1) + (1 << i2);
    for (int j = 0; j < idim; ++j) {
      int ibit = list1[j] & is;
      if (ibit != 0 && ibit != is) {
        int iexchg = list1[j] ^ is;
        int ia = iexchg & irght;
        int ib = (iexchg & ilft) / ihfbit;
        corr += x[j] * x[list2[ia].first + list2[ib].second];
      }
    }
    sxx[k] = corr / 4;
  }
}

void xcorr(int n, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2) {
  xcorr(n, npair, x.data(), sxx, list1, list2);
}

void xcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& sxx, std::vector<int>& list1,
           std::vector<std::pair<int, int>>& list2) {
  xcorr(n, npair, &x(0, xindex), sxx, list1, list2);
}

void zcorr(int n, std::vector<int> const& npair, const double *x,
           std::vector<double>& szz, std::vector<int>& list1) {
  int idim = list1.size();
  int nbond = npair.size() / 2;

  for (int k = 0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    if (i1 < 0 || i1 >= n || i2 < 0 || i2 >= n || i1 == i2) {
      std::cerr << " #(W02)# Wrong site number given to zcorr\n";
      return;
    }
    double corr = 0;
    int is = (1 << i1) + (1 << i2);
    for (int j = 0; j < idim; ++j) {
      int ibit = list1[j] & is;
      double factor;
      if (ibit == 0 || ibit == is) {
        factor = 1;
      } else {
        factor = -1;
      }
      corr += factor * x[j] * x[j];
    }
    szz[k] = corr / 4;
  }
}

void zcorr(int n, std::vector<int> const& npair, matrix_type const& x, int xindex,
           std::vector<double>& szz, std::vector<int>& list1) {
  zcorr(n, npair, &x(0, xindex), szz, list1);
}

void zcorr(int n, std::vector<int> const& npair, std::vector<double> const& x,
           std::vector<double>& szz, std::vector<int>& list1) {
  zcorr(n, npair, &x[0], szz, list1);
}

int orthg(matrix_type& ev, std::vector<double>& norm, int numvec) {
  if (numvec <= 1) {
    std::cerr << " #(W03)# Number of vectors is less than 2 in orthg\n";
    return -1;
  }
  if (norm.size() < numvec) norm.resize(numvec);
  int idim = ev.size1();
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
