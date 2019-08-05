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

#ifndef ROKKO_SCALAPACK_HPP
#define ROKKO_SCALAPACK_HPP

#include <rokko/cscalapack.h>

namespace rokko {
namespace scalapack {

inline int descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ictxt, int lld) {
  return cscalapack_descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld);
}
  
template<typename MATRIX, typename VECTOR>
int syev(char jobz, char uplo, MATRIX& a, int ia, int ja, VECTOR& w, MATRIX& z, int iz, int jz);

template<typename MATRIX, typename VECTOR>
int syev(char jobz, char uplo, MATRIX& a, int ia, int ja, VECTOR& w, MATRIX& z, int iz, int jz,
         VECTOR& work, int lwork);

} // end namespace scalapack
} // end namespace rokko

#include "scalapack/psyev.hpp"

#endif // ROKKO_SCALAPACK_HPP
