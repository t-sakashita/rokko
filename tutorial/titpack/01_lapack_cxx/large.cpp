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

#include "large.hpp"

int lnc1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
         std::vector<double> const& zrtio, int nvec, int iv, std::vector<double>& E,
         std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
         matrix_type& wk, std::vector<int> const& list1,
         std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  if (iv < 0 ||  iv >= idim) {
    std::cerr << " #(E06)# Incorrect iv given to lnc1\n";
    return -1;
  }
  if (nvec < 0 || nvec > 4) {
    std::cerr << " #(W06)# Wrong value given to nvec in lnc1\n"
              << "         Only the eigenvalues are calculated\n";
    nvec = 0;
  }
  return lnc1z(n, ipair, bondwt, zrtio, nvec, iv, E, alpha, beta, coeff, &wk(0,0), &wk(1,0),
               list1, list2);
}

int lnc1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, int nvec, int iv, std::vector<double>& E,
          std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
          double *v1, double *v0, std::vector<int> const& list1,
          std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  std::vector<int> iblock, isplit;
  matrix_type work(5, 150);
  double eps = 1e-10;
  int m, nsplit;
  double ebefor;

  // initialization
  alpha.clear();
  beta.clear();
  for (int i = 0; i < idim; ++i) {
    v0[i] = 0;
    v1[i] = 0;
  }
  v1[iv] = 1;
  
  datack(ipair, n);

  // alpha[0] and beta[0]
  double prdct = mltply(n, ipair, bondwt, zrtio, v1, v0, list1, list2);
  double alpha0 = prdct;
  alpha.push_back(alpha0);
  double beta0 = 0;
  for (int i = 0; i < idim; ++i) beta0 += (v0[i] - alpha0 * v1[i]) * (v0[i] - alpha0 * v1[i]);
  beta0 = std::sqrt(beta0);
  beta.push_back(beta0);
  
  // iteration  
  for (int i = 1; i < 150; ++i) {
    for (int j = 0; j < idim; ++j) {
      double temp1 = v1[j];
      double temp2 = (v0[j] - alpha0 * v1[j]) / beta0;
      v0[j] = -beta0 * temp1;
      v1[j] = temp2;
    }
    prdct = mltply(n, ipair, bondwt, zrtio, v1, v0, list1, list2);
    alpha0 = prdct;
    alpha.push_back(alpha0);
    beta0 = 0;
    for (int j = 0; j < idim; ++j) beta0 += (v0[j] - alpha0 * v1[j]) * (v0[j] - alpha0 * v1[j]);
    beta0 = std::sqrt(beta0);
    beta.push_back(beta0);
    if (beta[i] < 0.5e-30) {
      std::cerr << " #(E07)# Tridiagonalization unsuccessful in lnc1\n"
                << "         Beta(i) is too small at i=  " << i << std::endl;
      std::abort();
    }

    // convergence check
    if ((i + 1 > 20) && ((i + 1) % 5 == 0)) {
      boost::tie(m, nsplit) = bisec(alpha, beta, i + 1, E, 4, eps, iblock, isplit, &work(0,0));
      if (std::abs((ebefor - E[1]) / E[1]) < 1e-13) {
        if (nvec > 0) vec12(alpha, beta, i + 1, E, nvec, coeff, iblock, isplit, &work(0,0));
        return (i + 1);
      }
      ebefor = E[1];
    }
    if ((i + 1) == 20) {
      boost::tie(m, nsplit) = bisec(alpha, beta, 20, E, 4, eps, iblock, isplit, &work(0,0));
      ebefor = E[1];
    }
  }
  std::cerr << " #(W07)# lnc1 did not converge within 150 steps\n";
  return 150;
}

void lncv1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
           std::vector<double> const& zrtio, int nvec, int iv, std::vector<double> const& alpha,
           std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
           int itr, matrix_type& wk, std::vector<int> const& list1,
           std::vector<std::vector<int> > const& list2) {
  if (nvec <= 0 || nvec > 4) {
    std::cerr << " #(W08)# nvec given to lncv1 out of range\n";
    return;
  }
  int idim = list1.size();
  if (wk.size1() < 2 || wk.size2() != idim) wk.resize(2, idim);
  lncv1z(n, ipair, bondwt, zrtio, nvec, iv, alpha, beta, coeff, x, itr, &wk(0,0), &wk(1,0),
         list1, list2);
}

void lncv1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
            std::vector<double> const& zrtio, int nvec, int iv, std::vector<double> const& alpha,
            std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
            int itr, double *v1, double *v0, std::vector<int> const& list1,
            std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  if (x.size1() < nvec || x.size2() != idim) x.resize(nvec, idim);

  // initialization
  for (int i = 0; i < idim; ++i) {
    v0[i] = 0;
    v1[i] = 0;
  }
  v1[iv] = 1;
  
  for (int k = 0; k < nvec; ++k) {
    for (int i = 1; i < idim; ++i) x(k, i) = 0;
    x(k, iv) = coeff(k, 0);
  }

  // alpha(0) and beta(0)
  mltply(n, ipair, bondwt, zrtio, v1, v0, list1, list2);
  double alpha0 = alpha[0];
  double beta0 = beta[0];
  for (int k = 0; k < nvec; ++k)
    for (int j = 0; j < idim; ++j)
      x(k, j) += coeff(k, 1) * (v0[j] - alpha0 * v1[j]) / beta0;

  // iteration
  for (int i = 1; i < itr - 1; ++i) {
    for (int j = 0; j < idim; ++j) {
      double temp1 = v1[j];
      double temp2 = (v0[j] - alpha0 * v1[j]) / beta0;
      v0[j] = -beta0 * temp1;
      v1[j] = temp2;
    }
    mltply(n, ipair, bondwt, zrtio, v1, v0, list1, list2);
    alpha0 = alpha[i];
    beta0 = beta[i];
    for (int k = 0; k < nvec; ++k)
      for (int j = 1; j < idim; ++j)
        x(k, j) += coeff(k, i + 1) * (v0[j] - alpha0 * v1[j]) / beta0;
  }
  
  // normalization
  for (int k = 0; k < nvec; ++k) {
    double dnorm = 0;
    for (int j = 0; j < idim; ++j) dnorm += x(k, j) * x(k, j);
    dnorm = std::sqrt(dnorm);
    for (int j = 0; j < idim; ++j) x(k, j) /= dnorm;
  }
}

double mltply(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, double *v1, double *v0,
              std::vector<int> const& list1, std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  int ibond = ipair.size() / 2;

  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;

  double prdct = 0;
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[k * 2];
    int isite2 = ipair[k * 2 + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = 0.5 * bondwt[k] * zrtio[k];
    for (int j = 0; j < idim; ++j) {
      int ibit = list1[j] & is;
      double factor, offdg;
      if (ibit == 0 || ibit == is) {
        v0[j] -= wght * v1[j];
        factor = 1;
        offdg = 0;
      } else {
        v0[j] += wght * v1[j];
        int iexchg = list1[j] ^ is;
        int ia = iexchg & irght;
        int ib = (iexchg & ilft) / ihfbit;
        double temp = v1[list2[0][ia] + list2[1][ib]] * bondwt[k];
        v0[j] -= temp;
        factor = -1;
        offdg = -temp * v1[j];
      }
      prdct += - factor * wght * v1[j] * v1[j] + offdg;
    }
  }
  return prdct;
}

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, const double *x,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  int ibond = ipair.size() / 2;
  if (v.size1() < vindex || v.size2() != idim) v.resize(vindex, idim);

  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;
  
  double dnorm = 0;
  for (int i = 0; i < idim; ++i) dnorm += x[i] * x[i];
  if (dnorm < 1e-30) {
    std::cerr << " #(W09)# Null vector given to check1\n";
    return 0;
  }

  for (int i = 0; i < idim; ++i) v(vindex, i) = 0;
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[k * 2];
    int isite2 = ipair[k * 2 + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = bondwt[k];
    double diag = wght * 0.5 * zrtio[k];
    for (int i = 0; i < idim; ++i) {
      int ibit = list1[i] & is;
      if (ibit == 0 || ibit == is) {
        v(vindex, i) -= diag * x[i];
      } else {
        v(vindex, i) += diag * x[i];
        int iexchg = list1[i] ^ is;
        int ia = iexchg & irght;
        int ib = (iexchg & ilft) / ihfbit;
        v(vindex, i) -= x[list2[0][ia] + list2[1][ib]] * wght;
      }
    }
  }
  
  double prd = 0;
  for (int i = 0; i < idim; ++i) prd += v(vindex, i) * x[i];
  std::cout << "---------------------------- Information from check1\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v(vindex, i) / x[i];
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, matrix_type const& x, int xindex,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2) {
  return check1(n, ipair, bondwt, zrtio, &x(xindex, 0), v, vindex, list1, list2);
}

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, std::vector<double> const& x,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2) {
  return check1(n, ipair, bondwt, zrtio, &x[0], v, vindex, list1, list2);
}

void inv1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, double Eig, int iv, std::vector<double>& x,
          matrix_type& wk, std::vector<int> const& list1,
          std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  if (wk.size1() < 4 || wk.size2() != idim) wk.resize(4, idim);
  inv1z(n, ipair, bondwt, zrtio, Eig, iv, x, &wk(0,0), &wk(1,0), &wk(2,0), &wk(3,0), list1, list2);
}

void inv1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
           std::vector<double> const& zrtio, double Eig, int iv, std::vector<double>& x,
           double *b, double *p, double *r, double *y, std::vector<int> const& list1,
           std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  if (x.size() != idim) x.resize(idim);
  for (int i = 0; i < idim; ++i) b[i] = 0;
  b[iv]=1;

  for (int itr = 0; itr < 20; ++itr) {
    int iterat = cg1(n, ipair, bondwt, zrtio, Eig, x, b, p, r, y, list1, list2);
    if (iterat > idim) {
      std::cerr << " #(W10)# Iterat in cg1 exceeds idim or 500\n"
                << "         Approximate eigenvector returned\n"
                << "         Itration number in inv1 is " << itr << std::endl;
      return;
    }
    double xnorm = 0;
    for (int i =0; i < idim; ++i) xnorm += x[i] * x[i];
    xnorm = std::sqrt(xnorm);
    for (int i =0; i < idim; ++i) x[i] /= xnorm;
    double xb = 0;
    for (int i = 0; i < idim; ++i) xb += x[i] * b[i];
    if (std::abs(std::abs(xb) - 1) < 1e-12) {
      // std::cout << "       number of iterations in inv1 : " << (itr + 1) << std::endl;
      return;
    }
    for (int i = 0; i < idim; ++i) b[i] = x[i];
  }
  std::cerr << " #(W11)# inv1 did not converge\n";
  return;
}

int cg1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
        std::vector<double> const& zrtio, double Eig, std::vector<double>& x,
        double *b, double *p, double *r, double *y, std::vector<int> const& list1,
        std::vector<std::vector<int> > const& list2) {
  int idim = list1.size();
  int ibond = ipair.size() / 2;

  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;

  // initialization
  double bnorm = 0;
  for (int i=0; i < idim; ++i) {
    bnorm += b[i] * b[i];
    r[i] = b[i];
    p[i] = b[i];
    x[i] = 0;
  }

  // iteration
  for (int itr = 0; itr < std::min(500, idim); ++itr) {
    for (int i = 0; i < idim; ++i) y[i] = 0;
    for (int k = 0; k < ibond; ++k) {
      int isite1 = ipair[k * 2];
      int isite2 = ipair[k * 2 + 1];
      int is1 = 1 << isite1;
      int is2 = 1 << isite2;
      int is = is1 + is2;
      double eperbd = Eig / ibond;
      double wght = bondwt[k];
      double diag1 = wght * 0.5 * zrtio[k] + eperbd;
      double diag2 = -wght * 0.5 * zrtio[k] + eperbd;
      for (int i = 0; i < idim; ++i) {
        int ibit = list1[i] & is;
        if (ibit == 0 || ibit == is) {
          y[i] += (-diag1 * p[i]);
        } else {
          int iexchg = list1[i] ^ is;
          int ia = iexchg & irght;
          int ib = (iexchg & ilft) / ihfbit;
          y[i] += (-diag2 * p[i] - p[list2[0][ia] + list2[1][ib]] * wght);
        }
      }
    }
    double rp = 0;
    double yp = 0;
    for (int i = 0; i < idim; ++i) {
      rp += r[i] * p[i];
      yp += y[i] * p[i];
    }
    double alpha = rp / yp;
    double rnorm = 0;
    for (int i = 0; i < idim; ++i) {
      x[i] += alpha * p[i];
      rnorm += r[i] * r[i];
    }
    double rnorm2 = 0;
    for (int i = 0; i < idim; ++i) {
      r[i] -= alpha * y[i];
      rnorm2 += r[i] * r[i];
    }
    double beta = rnorm2 / rnorm;
    for (int i = 0; i < idim; ++i) p[i] = r[i] + beta * p[i];
    if ((itr + 1) % 5 == 0) {
      if (std::sqrt(rnorm2) < 1e-9 * std::sqrt(bnorm)) {
        // std::cout << "       number of iterations in cg1     : " << (itr + 1) << std::endl;
        return (itr + 1);
      }
    }
  }
  std::cerr << " #(Wxx)# cg1 did not converge\n";
  return std::min(500, idim);
}
