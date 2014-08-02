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

#ifndef TITPACK_LARGE_HPP
#define TITPACK_LARGE_HPP

#include "common.hpp"
#include <vector>

// The following variables are common to all routines
// n         @ lattice size
// ipair     @ pairs of sites connected by bonds
// bondwt    @ exchange interaction of each bond Jxy
// zrtio     @ ratio of Jz to Jxy

//
// eigenvalues by the Lanczos method
// --- dummy routine for simple working area
//

int lnc1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
         std::vector<double> const& zrtio, int nvec, int iv, std::vector<double>& E,
         std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
         matrix_type& wk, std::vector<int> const& list1,
         std::vector<std::vector<int> > const& list2);
// return valur # number of iterations required for convergence
// nvec      @ number of eigenvectors to calculate in lncvec
// iv        @ location of the nonzero element of the initial vector
// E         # eigenvalues
// wk        working areas
// list1,list2 @ spin configurations generated in 'sz'

//
// eigenvalues by the Lanczos method
//

int lnc1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, int nvec, int iv, std::vector<double>& E,
          std::vector<double>& alpha, std::vector<double>& beta, matrix_type& coeff,
          double *v1, double *v0, std::vector<int> const& list1,
          std::vector<std::vector<int> > const& list2);

//
// eigenvector by the Lanczos method
//
  
void lncv1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
           std::vector<double> const& zrtio, int nvec, int iv, std::vector<double> const& alpha,
           std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
           int itr, matrix_type& wk, std::vector<int> const& list1,
           std::vector<std::vector<int> > const& list2);
// nvec      @ number of eigenvectors to be calculated
// iv        @ location of the nonzero element of the initial vector
// x         # eigenvector
// itr     @ number of interations for convergence
// wk         working area
// list1,list2 @ spin configurations generated in 'sz'

//
// eigenvector by the Lanczos method
//

void lncv1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
            std::vector<double> const& zrtio, int nvec, int iv, std::vector<double> const& alpha,
            std::vector<double> const& beta, matrix_type const& coeff, matrix_type& x,
            int itr, double *v1, double *v0, std::vector<int> const& list1,
            std::vector<std::vector<int> > const& list2);

//
// matrix multiplication
//

double mltply(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, double *v1, double *v0,
              std::vector<int> const& list1, std::vector<std::vector<int> > const& list2);

//
// check of the eigenvector and eigenvalue
//
// x           @ eigenvector to be checked
// v           # H*x
// list1,list2 @ spin configurations generated in 'sz'

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, const double *x,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2);

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, matrix_type const& x, int xindex,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2);

double check1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
              std::vector<double> const& zrtio, std::vector<double> const& x,
              matrix_type& v, int vindex, std::vector<int> const& list1,
              std::vector<std::vector<int> > const& list2);

//
// inverse iteration
//   --- dummy routine for simple working area
//

void inv1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, double Eig, int iv, std::vector<double>& x,
          matrix_type& wk, std::vector<int> const& list1,
          std::vector<std::vector<int> > const& list2);

//
// inverse iteration
//
// b            working area for the rhs of (H-E(approx))*x=b
// p,r,y        working area used in the routine cg1

void inv1z(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
           std::vector<double> const& zrtio, double Eig, int iv, std::vector<double>& x,
           double *b, double *p, double *r, double *y, std::vector<int> const& list1,
           std::vector<std::vector<int> > const& list2);

//
// solution of linear equations -- cg method
//
// return value # number of iterations required for convergence
// Eig         @ eigenvalue
// x           # eigenvector
// b             working area for the rhs of (H-E(approx))*x=b
// p,r,y         working area used in the routine cg
// list1,list2 @ spin configurations generated in 'sz'

int cg1(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
        std::vector<double> const& zrtio, double Eig, std::vector<double>& x,
        double *b, double *p, double *r, double *y, std::vector<int> const& list1,
        std::vector<std::vector<int> > const& list2);

#endif
