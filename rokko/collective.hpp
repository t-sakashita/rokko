/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_COLLECTIVE_HPP
#define ROKKO_COLLECTIVE_HPP

#include <mpi.h>
#include "distributed_matrix.hpp"
#include "localized_matrix.hpp"
#include "pblas/pblas.h"

namespace rokko {

template<typename MATRIX_MAJOR>
void gather(rokko::distributed_matrix<MATRIX_MAJOR> const& from, double* to, int root) {
  if (!from.is_col_major()) {
    std::cerr << "Error (gather): matrix_row_major is not supported\n";
    MPI_Abort(MPI_COMM_WORLD,67);
    exit(67);
  }

  int ictxt;
  ROKKO_blacs_get(-1, 0, &ictxt);
  char char_grid_major = (from.get_grid().is_row_major() ? 'R' : 'C');
  ROKKO_blacs_gridinit(&ictxt, char_grid_major, from.get_nprow(), from.get_npcol());
  int m = from.get_m_global();
  int n = from.get_n_global();
  int rsrc = (from.get_grid().is_row_major() ? (root / from.get_npcol()) :
              (root % from.get_nprow()));
  int csrc = (from.get_grid().is_row_major() ? (root % from.get_npcol()) :
              (root / from.get_nprow()));
  int descFrom[9], descTo[9];
  int info;
  ROKKO_descinit(descFrom, m, n, from.get_mb(), from.get_nb(), 0, 0, ictxt, from.get_lld(), &info);
  ROKKO_descinit(descTo, m, n, m, n, rsrc, csrc, ictxt, m, &info);
  for (int j = 0; j < n; ++j)
    ROKKO_pdcopy(m, from.get_array_pointer(), 1, (j+1), descFrom, 1, to, 1, (j+1), descTo, 1);

  ROKKO_blacs_gridexit(&ictxt);
}

template<typename MATRIX_MAJOR>
void gather(rokko::distributed_matrix<MATRIX_MAJOR> const& from, localized_matrix<MATRIX_MAJOR>& to,
           int root) {
  gather(from, &to(0,0), root);
}

template<typename MATRIX_MAJOR>
void scatter(const double* from, distributed_matrix<MATRIX_MAJOR>& to, int root) {
  if (!to.is_col_major()) {
    std::cerr << "Error (scatter): matrix_row_major is not supported\n";
    MPI_Abort(MPI_COMM_WORLD,67);
    exit(67);
  }

  int ictxt;
  ROKKO_blacs_get(-1, 0, &ictxt);
  char char_grid_major = (to.get_grid().is_row_major() ? 'R' : 'C');
  ROKKO_blacs_gridinit(&ictxt, char_grid_major, to.get_nprow(), to.get_npcol());
  int m = to.get_m_global();
  int n = to.get_n_global();
  int rsrc = (to.get_grid().is_row_major() ? (root / to.get_npcol()) : (root % to.get_nprow()));
  int csrc = (to.get_grid().is_row_major() ? (root % to.get_npcol()) : (root / to.get_nprow()));
  int descFrom[9], descTo[9];
  int info;
  ROKKO_descinit(descFrom, m, n, m, n, rsrc, csrc, ictxt, m, &info);
  ROKKO_descinit(descTo, m, n, to.get_mb(), to.get_nb(), 0, 0, ictxt, to.get_lld(), &info);

  for (int j = 0; j < n; ++j)
    ROKKO_pdcopy(m, from, 1, (j+1), descFrom, 1, to.get_array_pointer(), 1, (j+1), descTo, 1);

  ROKKO_blacs_gridexit(&ictxt);
}

template<typename MATRIX_MAJOR>
void scatter(localized_matrix<MATRIX_MAJOR> const& from, distributed_matrix<MATRIX_MAJOR>& to,
            int root) {
  scatter(&from(0,0), to, root);
}

} // namespace rokko

#endif // ROKKO_COLLECTIVE_HPP
