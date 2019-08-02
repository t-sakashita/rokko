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

#ifndef ROKKO_BLACS_HPP
#define ROKKO_BLACS_HPP

#include <rokko/cblacs.h>

namespace rokko {
namespace blacs {

inline void barrier(int ictxt, char score) { cblacs_barrier(ictxt, score); }
  
inline int descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ictxt, int lld) {
  return cblacs_descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld);
}
  
inline void exit(int conti) { cblacs_exit(conti); }
  
inline void free_blacs_system_handle(int& ISysCxt) {  cblacs_free_blacs_system_handle(&ISysCxt); }

inline int get(int context, int request) { return cblacs_get(context, request); }
  
inline void gridexit(int& ictxt) { cblacs_gridexit(&ictxt); }

inline void gridinfo(int ictxt, int& nprow, int& npcol, int& myrow, int& mycol) {
  cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
}

inline void gridinit(int& ictxt, char order, int nprow, int npcol) {
  cblacs_gridinit(&ictxt, order, nprow, npcol);
}

inline void pinfo(int& mypnum, int& nprocs) { cblacs_pinfo(&mypnum, &nprocs); }
  
inline int sys2blacs_handle(const MPI_Comm& comm) {
  MPI_Fint comm_f = MPI_Comm_c2f(comm);
  return cblacs_sys2blacs_handle(&comm_f);
}

}
}

#endif // ROKKO_BLACS_HPP
