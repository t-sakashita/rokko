/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_FRANK_MATRIX_C_H
#define ROKKO_FRANK_MATRIX_C_H

#include <rokko/rokko.h>

#ifdef __cplusplus
extern "C"{
#endif

void rokko_frank_matrix_generate_distributed_matrix(rokko_distributed_matrix* matrix);

void rokko_frank_matrix_generate_localized_matrix(rokko_localized_matrix* matrix);

#ifdef __cplusplus
}
#endif

#endif
