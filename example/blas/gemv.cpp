/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <rokko/eigen3.hpp>
#include <rokko/blas.hpp>

int main() {
  std::size_t n = 4;
  Eigen::MatrixXd a(n, n);
  a << 0.959291425205444, 0.257508254123736, 0.243524968724989, 0.251083857976031,
    0.547215529963803, 0.840717255983663, 0.929263623187228, 0.616044676146639,
    0.138624442828679, 0.254282178971531, 0.349983765984809, 0.473288848902729,
    0.149294005559057, 0.814284826068816, 0.196595250431208, 0.351659507062997;
  Eigen::VectorXd x(n);
  x << 0.830828627896291, 0.585264091152724, 0.549723608291140, 0.917193663829810;
  Eigen::VectorXd y(n);
  y << 0.961898080855054, 0.00463422413406744, 0.774910464711502, 0.817303220653433;
  double alpha = 2.3;
  double beta = 0.5;

  std::cout << "a = [\n" << a << "\n]\n";
  std::cout << "x = [" << x.transpose() << "]^t\n";
  std::cout << "y = [" << y.transpose() << "]^t\n";

  rokko::blas::gemv(CblasNoTrans, alpha, a, x, 1, beta, y, 1);

  std::cout << alpha << " * a * x + " << beta << " * y = [" << y.transpose() << "]^t\n";
}

/* gemv.m

A = [0.959291425205444,0.257508254123736,0.243524968724989,0.251083857976031;0.547215529963803,0.840717255983663,0.929263623187228,0.616044676146639;0.138624442828679,0.254282178971531,0.349983765984809,0.473288848902729;0.149294005559057,0.814284826068816,0.196595250431208,0.351659507062997]
x = [0.830828627896291;0.585264091152724;0.549723608291140;0.917193663829810]
y = [0.961898080855054;0.00463422413406744;0.774910464711502;0.817303220653433]
2.3 * A * x + 0.5 * y

*/
