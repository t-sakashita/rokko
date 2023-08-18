/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2023 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <string>
#include <vector>


namespace rokko {

auto get_command_line_args(int argc, char** argv) {
  std::vector<std::string> names;

  for(auto num=1; num < argc; ++num) {
    names.emplace_back(argv[num]);
  }

  return names;
}

} // namespace rokko
