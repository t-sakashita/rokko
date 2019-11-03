/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_COLLECTIVE_HPP
#define ROKKO_COLLECTIVE_HPP

#include <mpi.h>
#include <stdexcept>
#include <rokko/blacs.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/cpblas.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>

namespace rokko {

template<typename T, typename MATRIX_MAJOR>
void gather(rokko::distributed_matrix<T, MATRIX_MAJOR> const& from, T* to, int root) {
  if (!from.is_col_major()) {
    throw std::invalid_argument("gather: matrix_row_major is not supported");
  }
  int ictxt = from.get_grid().get_blacs_context();
  int m = from.get_m_global();
  int n = from.get_n_global();
  int rsrc = (from.get_grid().is_row_major() ? (root / from.get_npcol()) :
              (root % from.get_nprow()));
  int csrc = (from.get_grid().is_row_major() ? (root % from.get_npcol()) :
              (root / from.get_nprow()));
  const int* descFrom = from.get_mapping().get_blacs_descriptor();
  int descTo[9];
  scalapack::descinit(descTo, m, n, m, n, rsrc, csrc, ictxt, m);
  for (int j = 0; j < n; ++j)
    cpblas_pdcopy(m, from.get_array_pointer(), 0, j, descFrom, 1, to, 0, j, descTo, 1);
}

template<typename T, typename MATRIX_MAJOR>
void gather(rokko::distributed_matrix<T, MATRIX_MAJOR> const& from,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>>& to, int root) {
  gather(from, &to(0,0), root);
}

template<typename T, typename MATRIX_MAJOR>
void scatter(const T* from, distributed_matrix<T, MATRIX_MAJOR>& to, int root) {
  if (!to.is_col_major()) {
    throw std::invalid_argument("scatter: matrix_row_major is not supported");
  }
  int ictxt = to.get_grid().get_blacs_context();
  int m = to.get_m_global();
  int n = to.get_n_global();
  int rsrc = (to.get_grid().is_row_major() ? (root / to.get_npcol()) : (root % to.get_nprow()));
  int csrc = (to.get_grid().is_row_major() ? (root % to.get_npcol()) : (root / to.get_nprow()));
  int descFrom[9];
  scalapack::descinit(descFrom, m, n, m, n, rsrc, csrc, ictxt, m);
  const int* descTo = to.get_mapping().get_blacs_descriptor();
  for (int j = 0; j < n; ++j)
    cpblas_pdcopy(m, from, 0, j, descFrom, 1, to.get_array_pointer(), 0, j, descTo, 1);
}

template<typename T, typename MATRIX_MAJOR>
void scatter(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> const& from,
  distributed_matrix<T, MATRIX_MAJOR>& to, int root) {
  scatter(&from(0,0), to, root);
}

} // namespace rokko

#endif // ROKKO_COLLECTIVE_HPP
