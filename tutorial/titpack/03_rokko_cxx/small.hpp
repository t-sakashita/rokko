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

#ifndef TITPACK_SMALL_HPP
#define TITPACK_SMALL_HPP

#include "common.hpp"
#include "hamiltonian.hpp"
#include <vector>
#include <tuple>

//
// matrix elements
//

void elm3(hamiltonian const& hop, Eigen::MatrixXd& elemnt) {
  elemnt.setZero();
  for (int k = 0; k < hop.num_bonds(); ++k) {
    int isite1, isite2;
    std::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = hop.bond_weight(k);
    double diag = 0.5 * wght * hop.z_ratio(k);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hop.dimension(); ++i) {
      int ibit = hop.config(i) & is;
      if (ibit == 0 || ibit == is) {
        elemnt(i,i) -= diag;
      } else {
        elemnt(i, i) += diag;
        int newcfg = hop.config2index(hop.config(i) ^ is);
        elemnt(i, newcfg) -= wght;
      }
    }
  }
}

template<typename T, typename MAJOR>
void elm3(hamiltonian const& hop, rokko::distributed_matrix<T, MAJOR>& elemnt) {
  elemnt.set_zeros();
  for (int k = 0; k < hop.num_bonds(); ++k) {
    int isite1, isite2;
    std::tie(isite1, isite2) = hop.site_pair(k);
    int is = (1 << isite1) + (1 << isite2);
    double wght = hop.bond_weight(k);
    double diag = 0.5 * wght * hop.z_ratio(k);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hop.dimension(); ++i) {
      if (elemnt.has_global_row_index(i)) {
        int ibit = hop.config(i) & is;
        if (ibit == 0 || ibit == is) {
          elemnt.update_global(i, i, -diag);
        } else {
          elemnt.update_global(i, i, +diag);
          int newcfg = hop.config2index(hop.config(i) ^ is);
          elemnt.update_global(i, newcfg, -wght);
        }
      }
    }
  }
}

void elm3_sx(subspace const& ss, int i1, int i2, Eigen::MatrixXd& elemnt) {
  elemnt.setZero();
  if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
    std::cerr << " #(W01)# Wrong site number given to xcorr\n";
    return;
  }
  int is = (1 << i1) + (1 << i2);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ss.dimension(); ++i) {
    int ibit = ss.config(i) & is;
    if (ibit != 0 && ibit != is) {
      int newcfg = ss.config2index(ss.config(i) ^ is);
      elemnt(i, newcfg) += 0.25;
    }
  }
}

template<typename T, typename MAJOR>
void elm3_sx(subspace const& ss, int i1, int i2, rokko::distributed_matrix<T, MAJOR>& elemnt) {
  elemnt.set_zeros();
  if (i1 < 0 || i1 >= ss.num_sites() || i2 < 0 || i2 >= ss.num_sites() || i1 == i2) {
    std::cerr << " #(W01)# Wrong site number given to xcorr\n";
    return;
  }
  int is = (1 << i1) + (1 << i2);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ss.dimension(); ++i) {
    if (elemnt.has_global_row_index(i)) {
      int ibit = ss.config(i) & is;
      if (ibit != 0 && ibit != is) {
        int newcfg = ss.config2index(ss.config(i) ^ is);
        elemnt.update_global(i, newcfg, 0.25);
      }
    }
  }
}

//
// check of the eigenvector and eigenvalue
//
// elemnt    @ nonzero elements
// x         @ eigenvector to be checked
// xindex
// return value: Hexpec <x*H*x>

template<typename MATRIX>
double check3(MATRIX const& elemnt, MATRIX const& x, int xindex, MATRIX& y) {
  int idim = elemnt.rows();
  double dnorm = 0;
  for (int j = 0; j < idim; ++j) {
    dnorm += x(j, xindex) * x(j, xindex);
  }
  if (dnorm < 1e-30) {
    std::cerr << " #(W18)# Null vector given to check3\n";
    return 0;
  }
  for (int i = 0; i < idim; ++i) {
    y(i, xindex) = 0;
    for (int j = 0; j < idim; ++j) y(i, xindex) += elemnt(i, j) * x(j, xindex);
  }
  double prd = 0;
  for (int i = 0; i < idim; ++i) prd += y(i, xindex) * x(i, xindex);
  std::cout << "---------------------------- Information from check3\n"
            << "<x*H*x> = "<< prd << std::endl
            << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (count % 4 == 0) std::cout << std::endl;
    std::cout << '\t' << y(i, xindex) / x(i, xindex);
  }
  std::cout << std::endl << "---------------------------------------------------\n";
  return prd;
}

template<typename MATRIX>
double check3_mpi(MATRIX const& elemnt, MATRIX const& x, int xindex, MATRIX& y) {
  int idim = elemnt.get_n_global();
  double dnorm = dot_product(x, false, xindex, x, false, xindex);
  product_v(1.0, elemnt, false, x, false, xindex, 0.0, y, false, xindex);
  double prd = dot_product(x, false, xindex, y, false, xindex) / dnorm;

  if (x.is_gindex({0, xindex})) {
    std::cout << "---------------------------- Information from check3\n"
              << "<x*H*x> = "<< prd << std::endl
              << "H*x(j)/x(j) (j=min(idim/3,13)-1,idim,max(1,idim/20))";
  }
  std::cout << std::flush;
  MPI_Barrier(elemnt.get_grid().get_comm());
  int count = 0;
  for (int i = std::min((int)(idim / 3), 13) - 1; i < idim; i += std::max(1,idim/20), ++count) {
    if (x.is_gindex({i, xindex})) {
      if (count % 4 == 0) std::cout << std::endl;
      std::cout << '\t' << y.get_global(i, xindex) / x.get_global(i, xindex);
    }
    std::cout << std::flush;
    MPI_Barrier(elemnt.get_grid().get_comm());
  }
  if (x.is_gindex({0, xindex})) {
    std::cout << std::endl << "---------------------------------------------------\n";
  }
  std::cout << std::flush;
  return prd;
}

template<typename MATRIX>
void xcorr3_mpi(subspace const& ss, std::vector<int> const& npair, MATRIX const& x, int xindex,
                std::vector<double>& sxx, MATRIX& sx, MATRIX& y) {
  double dnorm = dot_product(x, false, xindex, x, false, xindex);
  int nbond = npair.size() / 2;
  for (int k = 0; k < nbond; ++k) {
    int i1 = npair[k * 2];
    int i2 = npair[k * 2 + 1];
    elm3_sx(ss, i1, i2, sx);
    product_v(1.0, sx, false, x, false, xindex, 0.0, y, false, xindex);
    sxx[k] =dot_product(x, false, xindex, y, false, xindex) / dnorm;
  }
}

#endif
