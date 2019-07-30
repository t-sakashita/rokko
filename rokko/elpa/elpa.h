/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ELPA_H
#define ROKKO_ELPA_H

#ifdef __cplusplus
extern "C" {
#endif

#define complex   _Complex
#include <elpa/elpa.h>
#include <assert.h>

#define assert_elpa_ok(x) assert(x == ELPA_OK)

#ifdef __cplusplus
}
#endif

#endif // ROKKO_ELPA_H
