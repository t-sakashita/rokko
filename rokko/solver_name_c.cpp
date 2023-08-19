/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver_name.h>
#include <rokko/utility/solver_name.hpp>
#include <rokko/copy_string.hpp>

void rokko_split_solver_name(char* str, char** library_ptr, char** routine_ptr) {
  const auto [tmp_library, tmp_routine] = rokko::split_solver_name(str);
  *library_ptr = copy_string(tmp_library);
  *routine_ptr = copy_string(tmp_routine);
}
