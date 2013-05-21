// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_SELFADJOINTRANK2UPTADE_H
#define EIGEN_SELFADJOINTRANK2UPTADE_H

namespace Eigen { 

namespace internal {

/* Optimized selfadjoint matrix += alpha * uv' + conj(alpha)*vu'
 * It corresponds to the Level2 syr2 BLAS routine
 */

template<typename Scalar, typename Index, typename UType, typename VType, int UpLo>
struct selfadjoint_rank2_update_selector;

template<typename Scalar, typename Index, typename UType, typename VType>
struct selfadjoint_rank2_update_selector<Scalar,Index,UType,VType,Lower>
{
  static void run(Scalar* mat, Index stride, const UType& u, const VType& v, Scalar alpha)
  {
    const Index size = u.size();
    for (Index i=0; i<size; ++i)
    {
      Map<Matrix<Scalar,Dynamic,1> >(mat+stride*i+i, size-i) +=
                        (conj(alpha)  * conj(u.coeff(i))) * v.tail(size-i)
                      + (alpha * conj(v.coeff(i))) * u.tail(size-i);
    }
  }
};

template<typename Scalar, typename Index, typename UType, typename VType>
struct selfadjoint_rank2_update_selector<Scalar,Index,UType,VType,Upper>
{
  static void run(Scalar* mat, Index stride, const UType& u, const VType& v, Scalar alpha)
  {
    const Index size = u.size();
    for (Index i=0; i<size; ++i)
      Map<Matrix<Scalar,Dynamic,1> >(mat+stride*i, i+1) +=
                        (conj(alpha)  * conj(u.coeff(i))) * v.head(i+1)
                      + (alpha * conj(v.coeff(i))) * u.head(i+1);
  }
};

template<bool Cond, typename T> struct conj_expr_if
  : conditional<!Cond, const T&,
      CwiseUnaryOp<scalar_conjugate_op<typename traits<T>::Scalar>,T> > {};

} // end namespace internal

} // end namespace Eigen

#endif // EIGEN_SELFADJOINTRANK2UPTADE_H
