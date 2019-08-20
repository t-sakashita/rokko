/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGENEXA_HPP
#define ROKKO_EIGENEXA_HPP

#include <rokko/config.h>
#include <rokko/ceigenexa.h>

namespace rokko {
namespace eigenexa {

inline void init() { ceigenexa_init(); }
inline void init(MPI_Comm comm) { ceigenexa_init1(comm); }
inline void init(MPI_Comm comm, char grid_major) { ceigenexa_init2(comm, grid_major); }

inline void free() { ceigenexa_free(); }
inline void free(int flag) { ceigenexa_free1(flag); }

template<typename GRID>
std::pair<int, int> get_matdims(const GRID& grid, int n) {
  int nx, ny;
  ceigenexa_get_matdims(grid.get_nprow(), grid.get_npcol(), n, &nx, &ny);
  return std::make_pair(nx, ny);
}
  
template<typename MATRIX, typename VECTOR>
void eigen_s(MATRIX& a, VECTOR& w, MATRIX& z, int m_forward = 48, int m_backword = 128,
             char mode = 'A') {
  ceigenexa_eigen_s(a.get_m_global(), a.get_m_global(), a.get_array_pointer(), a.get_lld(),
                    storage(w), z.get_array_pointer(), z.get_lld(), m_forward, m_backword, mode);
}
  
template<typename MATRIX, typename VECTOR>
void eigen_s(MATRIX& a, VECTOR& w, int m_forward = 48, int m_backword = 128, char mode = 'N') {
  ceigenexa_eigen_s(a.get_m_global(), a.get_m_global(), a.get_array_pointer(), a.get_lld(),
                    storage(w), 0, 0, m_forward, m_backword, mode);
}
  
template<typename MATRIX, typename VECTOR>
void eigen_sx(MATRIX& a, VECTOR& w, MATRIX& z, int m_forward = 48, int m_backword = 128,
              char mode = 'A') {
  ceigenexa_eigen_sx(a.get_m_global(), a.get_m_global(), a.get_array_pointer(), a.get_lld(),
                     storage(w), z.get_array_pointer(), z.get_lld(), m_forward, m_backword, mode);
}
  
template<typename MATRIX, typename VECTOR>
void eigen_sx(MATRIX& a, VECTOR& w, int m_forward = 48, int m_backword = 128, char mode = 'N') {
  ceigenexa_eigen_sx(a.get_m_global(), a.get_m_global(), a.get_array_pointer(), a.get_lld(),
                     storage(w), 0, 0, m_forward, m_backword, mode);
}
  
} // end namespace eigenexa
} // end namespace rokko

#endif // ROKKO_EIGENEXA_HPP
