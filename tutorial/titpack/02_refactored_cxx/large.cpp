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

int lnc1(hamiltonian const& hop, int nvec, int iv, std::vector<double>& E,
         std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
         matrix_type& wk) {
  if (iv < 0 ||  iv >= hop.dimension()) {
    std::cerr << " #(E06)# Incorrect iv given to lnc1\n";
    return -1;
  }
  if (nvec < 0 || nvec > 4) {
    std::cerr << " #(W06)# Wrong value given to nvec in lnc1\n"
              << "         Only the eigenvalues are calculated\n";
    nvec = 0;
  }
  if (wk.size1() != hop.dimension() || wk.size1() < 2) wk.resize(hop.dimension(), 2);
  return lnc1z(hop, nvec, iv, E, alpha, beta, coeff, &wk(0,0), &wk(0,1));
}

int lnc1z(hamiltonian const& hop, int nvec, int iv, std::vector<double>& E,
          std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
          double *v1, double *v0) {
  std::vector<int> iblock, isplit;
  matrix_type work(5, 150);
  constexpr double eps = 1e-10;
  int m, nsplit;
  double ebefor;

  // initialization
  alpha.clear();
  beta.clear();
  for (int i = 0; i < hop.dimension(); ++i) {
    v0[i] = 0;
    v1[i] = 0;
  }
  v1[iv] = 1;
  
  // datack(ipair, n);

  // alpha[0] and beta[0]
  double prdct = hop.multiply(v1, v0);
  double alpha0 = prdct;
  alpha.emplace_back(alpha0);
  double beta0 = 0;
  for (int i = 0; i < hop.dimension(); ++i) beta0 += (v0[i] - alpha0 * v1[i]) * (v0[i] - alpha0 * v1[i]);
  beta0 = std::sqrt(beta0);
  beta.emplace_back(beta0);
  
  // iteration  
  for (int i = 1; i < 150; ++i) {
    for (int j = 0; j < hop.dimension(); ++j) {
      double temp1 = v1[j];
      double temp2 = (v0[j] - alpha0 * v1[j]) / beta0;
      v0[j] = -beta0 * temp1;
      v1[j] = temp2;
    }
    prdct = hop.multiply(v1, v0);
    alpha0 = prdct;
    alpha.emplace_back(alpha0);
    beta0 = 0;
    for (int j = 0; j < hop.dimension(); ++j) beta0 += (v0[j] - alpha0 * v1[j]) * (v0[j] - alpha0 * v1[j]);
    beta0 = std::sqrt(beta0);
    beta.emplace_back(beta0);
    if (beta[i] < 0.5e-30) {
      std::cerr << " #(E07)# Tridiagonalization unsuccessful in lnc1\n"
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
  std::cerr << " #(W07)# lnc1 did not converge within 150 steps\n";
  return 150;
}

void lncv1(hamiltonian const& hop, int nvec, int iv, std::vector<double> const& alpha,
           std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
           int itr, matrix_type& wk) {
  if (nvec <= 0 || nvec > 4) {
    std::cerr << " #(W08)# nvec given to lncv1 out of range\n";
    return;
  }
  if (wk.size1() != hop.dimension() || wk.size2() < 2) wk.resize(hop.dimension(), 2);
  lncv1z(hop, nvec, iv, alpha, beta, coeff, x, itr, &wk(0,0), &wk(0,1));
}

void lncv1z(hamiltonian const& hop, int nvec, int iv, std::vector<double> const& alpha,
            std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
            int itr, double *v1, double *v0) {
  if (x.size1() != hop.dimension() || x.size2() < nvec) x.resize(hop.dimension(), nvec);

  // initialization
  for (int i = 0; i < hop.dimension(); ++i) {
    v0[i] = 0;
    v1[i] = 0;
  }
  v1[iv] = 1;
  
  for (int k = 0; k < nvec; ++k) {
    for (int i = 0; i < hop.dimension(); ++i) x(i, k) = 0;
    x(iv, k) = coeff(0, k);
  }

  // alpha(0) and beta(0)
  hop.multiply(v1, v0);
  double alpha0 = alpha[0];
  double beta0 = beta[0];
  for (int k = 0; k < nvec; ++k)
    for (int j = 0; j < hop.dimension(); ++j)
      x(j, k) += coeff(1, k) * (v0[j] - alpha0 * v1[j]) / beta0;

  // iteration
  for (int i = 1; i < itr - 1; ++i) {
    for (int j = 0; j < hop.dimension(); ++j) {
      double temp1 = v1[j];
      double temp2 = (v0[j] - alpha0 * v1[j]) / beta0;
      v0[j] = -beta0 * temp1;
      v1[j] = temp2;
    }
    hop.multiply(v1, v0);
    alpha0 = alpha[i];
    beta0 = beta[i];
    for (int k = 0; k < nvec; ++k)
      for (int j = 1; j < hop.dimension(); ++j)
        x(j, k) += coeff(i + 1, k) * (v0[j] - alpha0 * v1[j]) / beta0;
  }
  
  // normalization
  for (int k = 0; k < nvec; ++k) {
    double dnorm = 0;
    for (int j = 0; j < hop.dimension(); ++j) dnorm += x(j, k) * x(j, k);
    dnorm = std::sqrt(dnorm);
    for (int j = 0; j < hop.dimension(); ++j) x(j, k) /= dnorm;
  }
}

double check1(hamiltonian const& hop, const double *x, matrix_type& v, int vindex) {
  int ibond = hop.num_bonds();
  if (v.size1() != hop.dimension() || v.size2() < vindex) v.resize(hop.dimension(), vindex);

  double dnorm = 0;
  for (int i = 0; i < hop.dimension(); ++i) dnorm += x[i] * x[i];
  if (dnorm < 1e-30) {
    std::cerr << " #(W09)# Null vector given to check1\n";
    return 0;
  }

  for (int i = 0; i < hop.dimension(); ++i) v(i, vindex) = 0;
  hop.multiply(x, &v(0, vindex));
  
  double prd = 0;
  for (int i = 0; i < hop.dimension(); ++i) prd += v(i, vindex) * x[i];
  std::cout << "---------------------------- Information from check1\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(hop.dimension() / 3), 13) - 1; i < hop.dimension();
       i += std::max(1,hop.dimension()/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << v(i, vindex) / x[i];
  }
  std::cout << std::endl
            << "---------------------------------------------------\n";
  return prd;
}

double check1(hamiltonian const& hop, matrix_type const& x, int xindex, matrix_type& v,
              int vindex) {
  return check1(hop, &x(0, xindex), v, vindex);
}

double check1(hamiltonian const& hop, std::vector<double> const& x, matrix_type& v, int vindex) {
  return check1(hop, x.data(), v, vindex);
}

void inv1(hamiltonian const& hop, double Eig, int iv, std::vector<double>& x, matrix_type& wk) {
  if (wk.size1() != hop.dimension() || wk.size2() < 4) wk.resize(hop.dimension(), 4);
  inv1z(hop, Eig, iv, x, &wk(0,0), &wk(0,1), &wk(0,2), &wk(0,3));
}

void inv1z(hamiltonian const& hop, double Eig, int iv, std::vector<double>& x, double *b,
           double *p, double *r, double *y) {
  if (x.size() != hop.dimension()) x.resize(hop.dimension());
  for (int i = 0; i < hop.dimension(); ++i) b[i] = 0;
  b[iv]=1;

  for (int itr = 0; itr < 20; ++itr) {
    int iterat = cg1(hop, Eig, x, b, p, r, y);
    if (iterat > hop.dimension()) {
      std::cerr << " #(W10)# Iterat in cg1 exceeds hop.dimension() or 500\n"
                << "         Approximate eigenvector returned\n"
                << "         Itration number in inv1 is " << itr << std::endl;
      return;
    }
    double xnorm = 0;
    for (int i =0; i < hop.dimension(); ++i) xnorm += x[i] * x[i];
    xnorm = std::sqrt(xnorm);
    for (int i =0; i < hop.dimension(); ++i) x[i] /= xnorm;
    double xb = 0;
    for (int i = 0; i < hop.dimension(); ++i) xb += x[i] * b[i];
    if (std::abs(std::abs(xb) - 1) < 1e-12) {
      // std::cout << "       number of iterations in inv1 : " << (itr + 1) << std::endl;
      return;
    }
    for (int i = 0; i < hop.dimension(); ++i) b[i] = x[i];
  }
  std::cerr << " #(W11)# inv1 did not converge\n";
  return;
}

int cg1(hamiltonian const& hop, double Eig, std::vector<double>& x, double *b, double *p,
        double *r, double *y) {
  // initialization
  double bnorm = 0;
  for (int i = 0; i < hop.dimension(); ++i) {
    bnorm += b[i] * b[i];
    r[i] = b[i];
    p[i] = b[i];
    x[i] = 0;
  }

  // iteration
  for (int itr = 0; itr < std::min(500, hop.dimension()); ++itr) {
    for (int i = 0; i < hop.dimension(); ++i) y[i] = 0;
    hop.multiply(p, y);
    for (int i = 0; i < hop.dimension(); ++i) y[i] = Eig * p[i] - y[i];
    double rp = 0;
    double yp = 0;
    for (int i = 0; i < hop.dimension(); ++i) {
      rp += r[i] * p[i];
      yp += y[i] * p[i];
    }
    double alpha = rp / yp;
    double rnorm = 0;
    for (int i = 0; i < hop.dimension(); ++i) {
      x[i] += alpha * p[i];
      rnorm += r[i] * r[i];
    }
    double rnorm2 = 0;
    for (int i = 0; i < hop.dimension(); ++i) {
      r[i] -= alpha * y[i];
      rnorm2 += r[i] * r[i];
    }
    double beta = rnorm2 / rnorm;
    for (int i = 0; i < hop.dimension(); ++i) p[i] = r[i] + beta * p[i];
    if ((itr + 1) % 5 == 0) {
      if (std::sqrt(rnorm2) < 1e-9 * std::sqrt(bnorm)) {
        // std::cout << "       number of iterations in cg1     : " << (itr + 1) << std::endl;
        return (itr + 1);
      }
    }
  }
  std::cerr << " #(Wxx)# cg1 did not converge\n";
  return std::min(500, hop.dimension());
}
