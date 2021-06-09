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

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#define complex   _Complex
//#include <elpa/elpa.h>
// **** begining of fixed elpa.h ****
#include <limits.h>
//#include <complex.h>

#include <elpa/elpa_version.h>

struct elpa_struct;
typedef struct elpa_struct *elpa_t;

struct elpa_autotune_struct;
typedef struct elpa_autotune_struct *elpa_autotune_t;


#include <elpa/elpa_constants.h>
#include <elpa/elpa_generated_c_api.h>
#include <elpa/elpa_generated.h>
//#include <elpa/elpa_generic.h>

const char *elpa_strerr(int elpa_error);
// **** end of fixed elpa.h ****
#include <assert.h>

#define assert_elpa_ok(x) assert(x == ELPA_OK)

#ifdef __cplusplus
}
#endif
