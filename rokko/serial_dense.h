/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SERIAL_DENSE_H
#define ROKKO_SERIAL_DENSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <rokko/localized_matrix.h>
#include <rokko/localized_vector.h>

enum {
  rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
};

#include <rokko/serial_dense_ev.h>


#ifdef __cplusplus
}
#endif

#endif /* ROKKO_SERIAL_DENSE_H */
