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

#include <rokko/ceigenexa.h>
#include <rokko/eigenexa/eigenexa_interface.h>

void ceigenexa_eigen_s(int n, int nvec, double *a, int lda, double *w, double *z, int ldz,
                       int m_forward, int m_backward, char mode) {
  EIGENEXA_eigen_s(&n, &nvec, a, &lda, w, z, &ldz, &m_forward, &m_backward, &mode);
}
