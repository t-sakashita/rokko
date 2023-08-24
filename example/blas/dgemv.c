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

#include <stdio.h>
#include <cblas.h>
#include <rokko/cmatrix.h>

int main() {
  int n = 4;
  double** a = alloc_dmatrix(n, n);
  mat_elem(a, 0, 0) = 0.959291425205444;
  mat_elem(a, 0, 1) = 0.257508254123736;
  mat_elem(a, 0, 2) = 0.243524968724989;
  mat_elem(a, 0, 3) = 0.251083857976031;
  mat_elem(a, 1, 0) = 0.547215529963803;
  mat_elem(a, 1, 1) = 0.840717255983663;
  mat_elem(a, 1, 2) = 0.929263623187228;
  mat_elem(a, 1, 3) = 0.616044676146639;
  mat_elem(a, 2, 0) = 0.138624442828679;
  mat_elem(a, 2, 1) = 0.254282178971531;
  mat_elem(a, 2, 2) = 0.349983765984809;
  mat_elem(a, 2, 3) = 0.473288848902729;
  mat_elem(a, 3, 0) = 0.149294005559057;
  mat_elem(a, 3, 1) = 0.814284826068816;
  mat_elem(a, 3, 2) = 0.196595250431208;
  mat_elem(a, 3, 3) = 0.351659507062997;
  double* x = alloc_dvector(n);
  x[0] = 0.830828627896291;
  x[1] = 0.585264091152724;
  x[2] = 0.549723608291140;
  x[3] = 0.917193663829810;
  double* y = alloc_dvector(n);
  y[0] = 0.961898080855054;
  y[1] = 0.00463422413406744;
  y[2] = 0.774910464711502;
  y[3] = 0.817303220653433;
  double alpha = 2.3;
  double beta = 0.5;

  printf("a: "); fprint_dmatrix(stdout, n, n, a);
  printf("x: "); fprint_dvector(stdout, n, x);
  printf("y: "); fprint_dvector(stdout, n, y);

  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, mat_ptr(a), n, vec_ptr(x), 1,
              beta, vec_ptr(y), 1);

  printf("%10.5f * a * x + %10.5f * y: ", alpha, beta); fprint_dvector(stdout, n, y);

  free_dmatrix(a);
  free_dvector(x);
  free_dvector(y);
}

/* gemv.m

A = [0.959291425205444,0.257508254123736,0.243524968724989,0.251083857976031;0.547215529963803,0.840717255983663,0.929263623187228,0.616044676146639;0.138624442828679,0.254282178971531,0.349983765984809,0.473288848902729;0.149294005559057,0.814284826068816,0.196595250431208,0.351659507062997]
x = [0.830828627896291;0.585264091152724;0.549723608291140;0.917193663829810]
y = [0.961898080855054;0.00463422413406744;0.774910464711502;0.817303220653433]
2.3 * A * x + 0.5 * y

*/
