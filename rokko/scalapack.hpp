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

#include <rokko/config.h>
#include <rokko/cscalapack.h>

namespace rokko {
namespace scalapack {

inline int descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ictxt, int lld) {
  return cscalapack_descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld);
}
  
template<typename MATRIX>
typename MATRIX::value_type pelget(char scope, char top, const MATRIX& A, int ia, int ja);

template<typename MATRIX>
void pelset(MATRIX& A, int ia, int ja, typename MATRIX::value_type alpha);

template<typename MATRIX>
typename MATRIX::value_type plange(char norm, const MATRIX& A);

template<typename MATRIX, typename VECTOR>
typename MATRIX::value_type plange(char norm, const MATRIX& A, VECTOR& work);

template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z);

template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z, VECTOR& work);

template<typename MATRIX, typename VECTOR>
int psyevd(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z);

template<typename MATRIX, typename VECTOR>
int psyevd(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z, VECTOR& work);

} // end namespace scalapack
} // end namespace rokko

#include "scalapack/pelget.hpp"
#include "scalapack/pelset.hpp"
#include "scalapack/plange.hpp"
#include "scalapack/psyev.hpp"
#include "scalapack/psyevd.hpp"

#endif // ROKKO_SCALAPACK_HPP
