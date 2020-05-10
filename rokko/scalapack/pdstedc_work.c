/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/cscalapack.h>
#include <rokko/scalapack/scalapack_interface.h>

int cscalapack_pdstedc_work(char compz, int n, double* d, double* e,
                            double* Q, int iq, int jq, const int* descQ,
                            double* work, int lwork, int* iwork, int liwork) {
  int iq_f = iq + 1;
  int jq_f = jq + 1;
  int info;
  SCALAPACK_pdstedc(&compz, &n, d, e, Q, &iq_f, &jq_f, descQ,
                    work, &lwork, iwork, &liwork, &info);
  return info;
}
