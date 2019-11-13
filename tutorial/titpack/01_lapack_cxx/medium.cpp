/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// C++ version of TITPACK Ver.2 by H. Nishimori

#include "medium.hpp"

void elm2(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, matrix_type& elemnt, i_matrix_type& loc,
          std::vector<int> const& list1, std::vector<std::pair<int, int> > const& list2) {
  int idim = list1.size();
  int ibond = ipair.size() / 2;
  int ic = ibond + 1;

  int ihf = (n + 1) / 2;
  int ihfbit = 1 << ihf;
  int irght = ihfbit - 1;
  int ilft = ((1 << n) - 1) ^ irght;

  datack(ipair, n);
  // initialization
  if (elemnt.size1() != idim || elemnt.size2() < ic) elemnt.resize(idim, ic);
  if (loc.size1() != idim || loc.size2() < ic) loc.resize(idim, ic);
  for (int i = 0; i < ic; ++i)
    for (int j = 0; j < idim; ++j)
      elemnt(j, i) = 0;

  // diagonal elements
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[k * 2];
    int isite2 = ipair[k * 2 + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = -bondwt[k] * zrtio[k] * 0.5;
    for (int i = 0; i < idim; ++i) {
      int ibit = list1[i] & is;
      if (ibit == 0 || ibit == is) {
        elemnt(i, ic - 1) += wght;
      } else {
        elemnt(i, ic - 1) -= wght;
      }
      loc(i, ic - 1) = i;
    }
  }
  
  // off-diagonal elements
  for (int k = 0; k < ibond; ++k) {
    int isite1 = ipair[k * 2];
    int isite2 = ipair[k * 2 + 1];
    int is1 = 1 << isite1;
    int is2 = 1 << isite2;
    int is = is1 + is2;
    double wght = -bondwt[k];
    for (int i = 0; i < idim; ++i) {
      int ibit = list1[i] & is;
      if (ibit == 0 || ibit == is) {
        elemnt(i, k) = 0;
        loc(i, k) = i;
      } else {
        elemnt(i, k) = wght;
        int iexchg = list1[i] ^ is;
        int ia = iexchg & irght;
        int ib = (iexchg & ilft) / ihfbit;
        loc(i, k) = list2[ia].first + list2[ib].second;
      }
    }
  }
}

int lnc2(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
         std::vector<double>& E, std::vector<double>& alpha, std::vector<double>& beta,
         matrix_type& coeff, matrix_type& wk) {
  int idim = elemnt.size1();
  if (iv < 0 ||  iv >= idim) {
    std::cerr << " #(E08)# Incorrect iv given to lnc2\n";
    return -1;
  }
  if (nvec < 0 || nvec > 4) {
    std::cerr << " #(W12)# Wrong value given to nvec in lnc2\n"
              << "         Only the eigenvalues are calculated\n";
    nvec = 0;
  }
  if (wk.size1() != idim || wk.size2() < 2) wk.resize(idim, 2);
  return lnc2z(elemnt, loc, nvec, iv, E, alpha, beta, coeff, &wk(0,0), &wk(0,1));
}

int lnc2z(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
          std::vector<double>& E, std::vector<double>& alpha, std::vector<double>& beta,
          matrix_type& coeff, double *v1, double *v0) {
  int idim = elemnt.size1();
  int ic = elemnt.size2();
  std::vector<int> iblock, isplit;
  matrix_type work(150, 5);
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
  
  // alpha[0] and beta[0]
  double prdct = 0;
  for (int k = 0; k < ic; ++k) {
    for (int j = 0; j < idim; ++j) {
      double temp = elemnt(j, k) * v1[loc(j, k)];
      v0[j] += temp;
      prdct += v1[j] * temp;
    }
  }
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
    double prdct = 0;
    for (int k = 0; k < ic; ++k) {
      for (int j = 0; j < idim; ++j) {
        double temp = elemnt(j, k) * v1[loc(j, k)];
        v0[j] += temp;
        prdct += v1[j] * temp;
      }
    }
    alpha0 = prdct;
    alpha.push_back(alpha0);
    beta0 = 0;
    for (int j = 0; j < idim; ++j) beta0 += (v0[j] - alpha0 * v1[j]) * (v0[j] - alpha0 * v1[j]);
    beta0 = std::sqrt(beta0);
    beta.push_back(beta0);
    if (beta[i] < 0.5e-30) {
      std::cerr << " #(E09)# Tridiagonalization unsuccessful in lnc2\n"
                << "         Beta(i) is too small at i=  " << i << std::endl;
      std::abort();
    }

    // convergence check
    if ((i + 1 > 20) && ((i + 1) % 5 == 0)) {
      std::tie(m, nsplit) = bisec(alpha, beta, i + 1, E, 4, eps, iblock, isplit, &work(0,0));
      if (std::abs((ebefor - E[1]) / E[1]) < 1e-13) {
        if (nvec > 0) vec12(alpha, beta, i + 1, E, nvec, coeff, iblock, isplit, &work(0,0));
        return (i + 1);
      }
      ebefor = E[1];
    }
    if ((i + 1) == 20) {
      std::tie(m, nsplit) = bisec(alpha, beta, 20, E, 4, eps, iblock, isplit, &work(0,0));
      ebefor = E[1];
    }
  }
  std::cerr << " #(W13)# lnc2 did not converge within 150 steps\n";
  return 150;
}

void lncv2(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
           std::vector<double> const& alpha, std::vector<double> const& beta,
           matrix_type const& coeff, matrix_type& x, int itr, matrix_type& wk) {
  if (nvec <= 0 || nvec > 4) {
    std::cerr << " #(W14)# nvec given to lncv2 out of range\n";
    return;
  }
  int idim = elemnt.size1();
  if (wk.size1() != idim || wk.size2() < 2) wk.resize(idim, 2);
  lncv2z(elemnt, loc, nvec, iv, alpha, beta, coeff, x, itr, &wk(0,0), &wk(0,1));
}

void lncv2z(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
            std::vector<double> const& alpha, std::vector<double> const& beta,
            matrix_type const& coeff, matrix_type& x, int itr, double *v1, double *v0) {
  int idim = elemnt.size1();
  int ic = elemnt.size2();
  if (x.size1() != idim || x.size2() < nvec) x.resize(idim, nvec);

  // initialization
  for (int i = 0; i < idim; ++i) {
    v0[i] = 0;
    v1[i] = 0;
  }
  v1[iv] = 1;
  
  for (int k = 0; k < nvec; ++k) {
    for (int i = 1; i < idim; ++i) x(i, k) = 0;
    x(iv, k) = coeff(0, k);
  }

  // alpha(0) and beta(0)
  double prdct = 0;
  for (int k = 0; k < ic; ++k) {
    for (int j = 0; j < idim; ++j) {
      double temp = elemnt(j, k) * v1[loc(j, k)];
      v0[j] += temp;
      prdct += v1[j] * temp;
    }
  }
  double alpha0 = alpha[0];
  double beta0 = beta[0];
  for (int k = 0; k < nvec; ++k)
    for (int j = 0; j < idim; ++j)
      x(j, k) += coeff(1, k) * (v0[j] - alpha0 * v1[j]) / beta0;

  // iteration
  for (int i = 1; i < itr - 1; ++i) {
    for (int j = 0; j < idim; ++j) {
      double temp1 = v1[j];
      double temp2 = (v0[j] - alpha0 * v1[j]) / beta0;
      v0[j] = -beta0 * temp1;
      v1[j] = temp2;
    }
    double prdct = 0;
    for (int k = 0; k < ic; ++k) {
      for (int j = 0; j < idim; ++j) {
        double temp = elemnt(j, k) * v1[loc(j, k)];
        v0[j] += temp;
        prdct += v1[j] * temp;
      }
    }
    alpha0 = alpha[i];
    beta0 = beta[i];
    for (int k = 0; k < nvec; ++k)
      for (int j = 1; j < idim; ++j)
        x(j, k) += coeff(i + 1, k) * (v0[j] - alpha0 * v1[j]) / beta0;
  }
  
  // normalization
  for (int k = 0; k < nvec; ++k) {
    double dnorm = 0;
    for (int j = 0; j < idim; ++j) dnorm += x(j, k) * x(j, k);
    dnorm = std::sqrt(dnorm);
    for (int j = 0; j < idim; ++j) x(j, k) /= dnorm;
  }
}

double check2(matrix_type const& elemnt, i_matrix_type const& loc, const double *x,
              matrix_type& v, int vindex) {
  int idim = elemnt.size1();
  int ic = elemnt.size2();
  if (v.size1() != idim || v.size2() < vindex) v.resize(idim, vindex);

  double dnorm = 0;
  for (int i = 0; i < idim; ++i) dnorm += x[i] * x[i];
  if (dnorm < 1e-30) {
    std::cerr << " #(W15)# Null vector given to check2\n";
    return 0;
  }

  for (int i = 0; i < idim; ++i) v(i, vindex) = 0;
  for (int k = 0; k < ic; ++k)
    for (int j = 0; j < idim; ++j)
      v(j, vindex) += elemnt(j, k)  * x[loc(j, k)];
  
  double prd = 0;
  for (int i = 0; i < idim; ++i) prd += v(i, vindex) * x[i];
  std::cout << "---------------------------- Information from check2\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v(i, vindex) / x[i];
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}

double check2(matrix_type const& elemnt, i_matrix_type const& loc, matrix_type const& x, int xindex,
              matrix_type& v, int vindex) {
  return check2(elemnt, loc, &x(0, xindex), v, vindex);
}

double check2(matrix_type const& elemnt, i_matrix_type const& loc, std::vector<double> const& x,
              matrix_type& v, int vindex) {
  return check2(elemnt, loc, &x[0], v, vindex);
}

void inv2(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, int iv,
          std::vector<double>& x, matrix_type& wk) {
  int idim = elemnt.size1();
  if (wk.size1() != idim || wk.size2() < 4) wk.resize(idim, 4);
  inv2z(elemnt, loc, Eig, iv, x, &wk(0,0), &wk(0,1), &wk(0,2), &wk(0,3));
}

void inv2z(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, int iv,
           std::vector<double>& x, double *b, double *p, double *r, double *y) {
  int idim = elemnt.size1();
  int ic = elemnt.size2();
  if (x.size() != idim) x.resize(idim);
  for (int i = 0; i < idim; ++i) b[i] = 0;
  b[iv]=1;

  for (int itr = 0; itr < 20; ++itr) {
    int iterat = cg2(elemnt, loc, Eig, x, b, p, r, y);
    if (iterat > idim) {
      std::cerr << " #(W16)# Iterat in cg2 exceeds idim or 500\n"
                << "         Approximate eigenvector returned\n"
                << "         Itration number in inv2 is " << itr << std::endl;
      return;
    }
    double xnorm = 0;
    for (int i =0; i < idim; ++i) xnorm += x[i] * x[i];
    xnorm = std::sqrt(xnorm);
    for (int i =0; i < idim; ++i) x[i] /= xnorm;
    double xb = 0;
    for (int i = 0; i < idim; ++i) xb += x[i] * b[i];
    if (std::abs(std::abs(xb) - 1) < 1e-12) {
      // std::cout << "       number of iterations in inv2 : " << (itr + 1) << std::endl;
      return;
    }
    for (int i = 0; i < idim; ++i) b[i] = x[i];
  }
  std::cerr << " #(W11)# inv2 did not converge\n";
  return;
}

int cg2(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, std::vector<double>& x,
        double *b, double *p, double *r, double *y) {
  int idim = elemnt.size1();
  int ic = elemnt.size2();

  // initialization
  double bnorm = 0;
  for (int i = 0; i < idim; ++i) {
    bnorm += b[i] * b[i];
    r[i] = b[i];
    p[i] = b[i];
    x[i] = 0;
  }

  // iteration
  for (int itr = 0; itr < std::min(500, idim); ++itr) {
    for (int i = 0; i < idim; ++i) y[i] = -Eig * p[i];
    for (int k = 0; k < ic; ++k)
      for (int i = 0; i < idim; ++i)
        y[i] += elemnt(i, k) * p[loc(i, k)];
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
        // std::cout << "       number of iterations in cg2     : " << (itr + 1) << std::endl;
        return (itr + 1);
      }
    }
  }
  std::cerr << " #(Wxx)# cg2 did not converge\n";
  return std::min(500, idim);
}
