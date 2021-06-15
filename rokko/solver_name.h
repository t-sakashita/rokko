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

#pragma once

#ifdef __cplusplus
extern "C" {
#endif


//void rokko_split_solver_name(char* str, char* library, char* routine);

void rokko_split_solver_name(char* str, char** library_ptr, char** routine_ptr);

#ifdef __cplusplus
}
#endif
