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

#ifndef ROKKO_BLACS_UTILITY_ROUTINES_HPP
#define ROKKO_BLACS_UTILITY_ROUTINES_HPP

#include <rokko/Cblacs.h>
#include <rokko/scalapack.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {
namespace blacs {

template<typename MATRIX_MAJOR>
char set_grid(int& ictxt, distributed_matrix<double, MATRIX_MAJOR>& mat) {
  char char_grid_major;
  if(mat.get_grid().is_row_major())  char_grid_major = 'R';
  else char_grid_major = 'C';
  gridinit(ictxt, char_grid_major, mat.get_grid().get_nprow(), mat.get_grid().get_npcol());
  return char_grid_major;
}

template<typename MATRIX_MAJOR>
void set_desc(int ictxt, distributed_matrix<double, MATRIX_MAJOR>& mat, int desc[9]) {
  int lld = (mat.get_lld() != 0) ? mat.get_lld() : 1;
  int info = scalapack::descinit(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(),
                                 mat.get_nb(), 0, 0, ictxt, lld);
  int nprow, npcol, myrow, mycol;
  if (info) {
    std::cerr << "error info=" << info << " at descinit function of descA " << "mA="
              << mat.get_m_local() << "  nA=" << mat.get_n_local() << "  lld=" << mat.get_lld()
              << "." << std::endl;
    gridinfo(ictxt, nprow, npcol, myrow, mycol);
    std::cerr << "cblacs_gridinfo: " << "nprow=" << nprow << " npcol=" << npcol
	      << " myrow=" << myrow << " mycol=" << mycol << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }
}

} // namespace blacs
} // namespace rokko

#endif // ROKKO_BLACS_UTILITY_ROUTINES_HPP

