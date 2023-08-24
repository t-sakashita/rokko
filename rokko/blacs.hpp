/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <vector>
#include <rokko/Cblacs.h>

namespace rokko {
namespace blacs {

inline void barrier(int ictxt, char score) { Cblacs_barrier(ictxt, &score); }

inline MPI_Comm blacs2sys_handle(int BlacsCtxt) { return Cblacs2sys_handle(BlacsCtxt); }

inline void exit(int NotDone) { Cblacs_exit(NotDone); }

inline void free_blacs_system_handle(int ISysCxt) { Cfree_blacs_system_handle(ISysCxt); }

inline void get(int ConTxt, int what, int& val) { Cblacs_get(ConTxt, what, &val); }

inline void gridexit(int ConTxt) { Cblacs_gridexit(ConTxt); }

inline void gridinfo(int ConTxt, int& nprow, int& npcol, int& myrow, int& mycol) {
  Cblacs_gridinfo(ConTxt, &nprow, &npcol, &myrow, &mycol);
}

inline void gridinit(int& ConTxt, char order, int nprow, int npcol) {
  Cblacs_gridinit(&ConTxt, &order, nprow, npcol);
}

inline void gridmap(int& ConTxt, const std::vector<int>& usermap, int ldup,
                    int nprow0, int npcol0) {
  Cblacs_gridmap(&ConTxt, const_cast<int*>(usermap.data()), ldup, nprow0, npcol0);
}

inline int sys2blacs_handle(const MPI_Comm& comm) {
  return Csys2blacs_handle(comm);
}

} // end namespace blacs
} // end namespace rokko
