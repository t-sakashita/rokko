/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/xblas.hpp>
#include <rokko/lapackx.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>
#include <boost/lexical_cast.hpp>

int imin(int x, int y) { return (x < y) ? x : y; }
int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int m = 3;
  int n = 5;
  if (argc > 2) {
    m = boost::lexical_cast<int>(argv[1]);
    n = boost::lexical_cast<int>(argv[2]);
  }
  int r = imin(m, n);
  
  // generate matrix
  rokko::dlmatrix a(m, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      a(i, j) = n - imax(i, j);
    }
  }
  std::cout << "Matrix A: " << rows(a) << ' ' << cols(a) << std::endl
            << a << std::endl;

  // singular value decomposition
  rokko::dlvector s(r);
  rokko::dlmatrix u(m, r), vt(r, n);
  {
    rokko::dlmatrix t(a);
    rokko::dlvector work(imax(3 * r + imax(m, n), 5 * r));
    int info = rokko::lapackx::gesvd_work('S', 'S', t, s, u, vt, work);
  }
  std::cout << "Matrix U: " << rows(u) << ' ' << cols(u) << std::endl
            << u << std::endl;
  std::cout << "Singular values: " << size(s) << std::endl
            << s << std::endl;
  std::cout << "Matrix Vt: " << rows(vt) << ' ' << cols(vt) << std::endl
            << vt << std::endl;

  // orthogonality check
  {
    rokko::dlmatrix t(r, r);
    rokko::xblas::gemm(CblasTrans, CblasNoTrans, 1, u, u, 0, t);
    for (int i = 0; i < r; ++i) t(i, i) -= 1;
    double norm2 = 0;
    for (int j = 0; j < r; ++j) {
      for (int i = 0; i < r; ++i) {
        norm2 = t(i, j) * t(i, j);
      }
    }
    std::cout << "|| U^t U - I ||^2 = " << norm2 << std::endl;
    if (norm2 > 1e-16) {
      std::cerr << "Error: orthogonality check" << std::endl;
      exit(255);
    }
  }
  {
    rokko::dlmatrix t(r, r);
    rokko::xblas::gemm(CblasNoTrans, CblasTrans, 1, vt, vt, 0, t);
    for (int i = 0; i < r; ++i) t(i, i) -= 1;
    double norm2 = 0;
    for (int j = 0; j < r; ++j) {
      for (int i = 0; i < r; ++i) {
        norm2 = t(i, j) * t(i, j);
      }
    }
    std::cout << "|| V V^t - I ||^2 = " << norm2 << std::endl;
    if (norm2 > 1e-16) {
      std::cerr << "Error: orthogonality check" << std::endl;
      exit(255);
    }
  }

  // solution check
  {
    rokko::dlmatrix t(m, r);
    rokko::xblas::gemm(CblasNoTrans, CblasTrans, 1, a, vt, 0, t);
    rokko::dlmatrix w(r, r);
    rokko::xblas::gemm(CblasTrans, CblasNoTrans, 1, u, t, 0, w);
    for (int i = 0; i < r; ++i) w(i, i) -= s(i);
    double norm2 = 0;
    for (int j = 0; j < r; ++j) {
      for (int i = 0; i < r; ++i) {
        norm2 = w(i, j) * w(i, j);
      }
    }
    std::cout << "|| U^t A V - diag(S) ||^2 = " << norm2 << std::endl;
    if (norm2 > 1e-16) {
      std::cerr << "Error: solution check" << std::endl;
      exit(255);
    }
  }

  return 0;
}
