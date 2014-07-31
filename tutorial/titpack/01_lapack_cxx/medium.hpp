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

#ifndef TITPACK_MEDIUM_HPP
#define TITPACK_MEDIUM_HPP

#include "common.hpp"
#include <vector>

// The following variables are common to all routines
// elemnt    @ nonzero elements
// loc       @ location of nonzero elements

//
// matrix elements for general bond weights
//
// n           @ lattice size
// ipair       @ pairs of sites connected by bonds
// bondwt      @ exchange interaction of each bond Jxy
// zrtio       @ ratio of Jz to Jxy
// list1,list2 @ @spin configurations generated in 'sz'

void elm2(int n, std::vector<int> const& ipair, std::vector<double> const& bondwt,
          std::vector<double> const& zrtio, matrix_type& elemnt, i_matrix_type& loc,
          std::vector<int> const& list1, std::vector<std::vector<int> > const& list2);

//
// eigenvalues by the Lanczos method
// --- dummy routine for simple working area
//

int lnc2(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
         std::vector<double>& E, std::vector<double>& alpha, std::vector<double>& beta,
         matrix_type& coeff, matrix_type& wk);
// return value # number of iterations required for convergence
// nvec      @ number of eigenvectors to calculate in lncvec
// iv        @ location of the nonzero element of the initial vector
// E         # eigenvalues
// wk        working areas

//
// eigenvalues by the Lanczos method
//

int lnc2z(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
          std::vector<double>& E, std::vector<double>& alpha, std::vector<double>& beta,
          matrix_type& coeff, double *v1, double *v0);

//
// eigenvector by the Lanczos method
//
  
void lncv2(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
           std::vector<double> const& alpha, std::vector<double> const& beta,
           matrix_type const& coeff, matrix_type& x, int itr, matrix_type& wk);
// nvec      @ number of eigenvectors to be calculated
// iv        @ location of the nonzero element of the initial vector
// x         # eigenvector
// itr     @ number of interations for convergence
// wk         working area

//
// eigenvector by the Lanczos method
//

void lncv2z(matrix_type const& elemnt, i_matrix_type const& loc, int nvec, int iv,
            std::vector<double> const& alpha, std::vector<double> const& beta,
            matrix_type const& coeff, matrix_type& x, int itr, double *v1, double *v0);

//
// check of the eigenvector and eigenvalue
//
// x           @ eigenvector to be checked
// v           # H*x

double check2(matrix_type const& elemnt, i_matrix_type const& loc, const double *x,
              matrix_type& v, int vindex);

double check2(matrix_type const& elemnt, i_matrix_type const& loc, matrix_type const& x, int xindex,
              matrix_type& v, int vindex);

double check2(matrix_type const& elemnt, i_matrix_type const& loc, std::vector<double> const& x,
              matrix_type& v, int vindex);

//
// inverse iteration
//   --- dummy routine for simple working area
//

void inv2(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, int iv,
          std::vector<double>& x, matrix_type& wk);

//
// inverse iteration
//
// b            working area for the rhs of (H-E(approx))*x=b
// p,r,y        working area used in the routine cg1

void inv2z(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, int iv,
           std::vector<double>& x, double *b, double *p, double *r, double *y);

//
// solution of linear equations -- cg method
//
// return value # number of iterations required for convergence
// Eig         @ eigenvalue
// x           # eigenvector
// b             working area for the rhs of (H-E(approx))*x=b
// p,r,y         working area used in the routine cg
// list1,list2 @ spin configurations generated in 'sz'

int cg2(matrix_type const& elemnt, i_matrix_type const& loc, double Eig, std::vector<double>& x,
        double *b, double *p, double *r, double *y);

#endif
